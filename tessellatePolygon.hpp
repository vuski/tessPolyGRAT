#ifndef VV_TESSELLATE_POLYGON_HPP
#define VV_TESSELLATE_POLYGON_HPP

#include "spdlog/spdlog.h"
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>



using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::unordered_map;




namespace vv {

class TessPolyGRAT {

public:



	void init(double gridUnitGlobal_, uint32_t basex = 740000,
		uint32_t maxx = 1390000,
		uint32_t basey = 1450000,
		uint32_t maxy = 2070000,
		bool isDebugMode = false)
	{
		gridUnitGlobal = gridUnitGlobal_;
		BASEX = basex;
		MAXX = maxx;
		BASEY = basey;
		MAXY = maxy;
		debug01 = isDebugMode;
	}


	//�������� �߰��� �������� �׼����̼� �Ѵ�.
	template <typename VECTOR2D>
	bool tessellatePolygon(vector<vector<VECTOR2D>>& in_polygons,
		unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& out_tessellated,
		vector<VECTOR2D>& out_vertex)
	{
		uint32_t gridunit = (uint32_t)round(gridUnitGlobal);

		//�׸���� �����ϴ� ����
		vector<IntersectPt> axis_x;

		//Ȧ�� ������ �̱� �����￡ �׸������ ���������� �߰��Ѵ�.
		bool isOK = addIntersectedPointsToPolygonByTriangleGrid(in_polygons, axis_x, gridunit);
		if (!isOK) return false;

		//���������� �߰��� �������� �����Ѵ�. 
		//����� ��������� ���Ͱ� ��ȯ�Ǿ�� �Ѵ�. ������� ����
		//��ġ�� ������ �����Ƿ� ���ؽ��� �ε����� �����ؾ� �Ѵ�.

		//���� �������� ���ؽ��� �δ��̷�Ʈ�� ��ȯ
		//key�� �׸��� ���ڵ� ��ȣ�� �ﰢ�� �ڸ� ��
		unordered_map<uint32_t, vector<GRID_STRG>> gridStorage;
		unordered_map<uint32_t, uint32_t> gridPointMap; //��ǥ xy�ε����� ������ vertex vector�� �ε����� ��ȯ�Ѵ�.

		uint32_t xrange = (MAXX - BASEX) / gridunit;
		uint32_t yrange = (MAXY - BASEY) / gridunit;


		//�ﰢ�� �׸��� ������, �������� ���׸�Ʈ�� �����ؼ� �ִ´�.
		distributeSegmentsToGrid(in_polygons, gridStorage, out_vertex, gridPointMap, gridunit, xrange, axis_x);
		//SPDLOG_DEBUG("���׸�Ʈ ��� �Ϸ�");


		//������� ����, �׸��庰 ��������� ������ ��� ����
		//���� ���ҵ� �׸��忡 ����� segment ���� �̾��ش�.
		//SPDLOG_DEBUG("��������� ���� �Ҵ� ��. ���� ����");
		constructPolygonFromSegments(out_tessellated, gridStorage, out_vertex, gridPointMap, gridunit, xrange);
		//SPDLOG_DEBUG("������ �籸�� �Ϸ�");


		//���������� �߰��� �� �������� ä���.
		fillTrianglesInside(out_tessellated, gridStorage, axis_x, out_vertex, gridPointMap, gridunit, xrange);



		return true;
	}


	template <typename VECTOR2D>
	void writeResult(unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& tessellated,
		vector<VECTOR2D>& vertex, string fileName, uint32_t gridsize)
	{
		uint32_t xrange = (MAXX - BASEX) / gridsize;

		auto startTime = std::chrono::high_resolution_clock::now();
		SPDLOG_DEBUG("tessellated geojson ���� ���� ����.....");
		std::stringstream header(std::stringstream::out | std::stringstream::binary);
		header << "{\"type\": \"FeatureCollection\"," << endl;
		header << "\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:EPSG::5179\" } }," << endl;
		header << "\"features\":[" << endl;
		string headerStr = header.str();

		auto myfile = std::fstream(fileName, std::ios::out | std::ios::binary);
		myfile.write(header.str().c_str(), header.str().length());

		std::stringstream str(std::stringstream::out | std::stringstream::binary);
		str.precision(20);

		unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>::iterator iter;
		for (iter = tessellated.begin(); iter != tessellated.end(); )
		{
			uint32_t gridxy = iter->first;
			uint32_t gridx = (gridxy % xrange) * gridsize + BASEX;
			uint32_t gridy = (gridxy / xrange) * gridsize + BASEY;
			vector<vector<vector<uint32_t>>>& multiPolygons = iter->second;
			bool isAvailable = false;
			str.str("");

			for (uint32_t i = 0; i < multiPolygons.size(); i++)
			{
				vector<vector<uint32_t>>& polygons = multiPolygons[i];

				str << "{ \"type\": \"Feature\", \"properties\": {";
				str << "\"gridx\":" << gridx << ", \"gridy\":" << gridy;
				str << " }, \"geometry\": { \"type\": \"Polygon\", \"coordinates\": [ ";


				for (uint32_t j = 0; j < polygons.size(); j++)
				{
					vector<uint32_t>& polygon = polygons[j];
					str << "[";
					for (uint32_t k = 0; k < polygon.size(); k++)
					{
						uint32_t& pt = polygon[k];
						str << "[" << vertex[pt].x << "," << vertex[pt].y;
						if (k == polygon.size() - 1) str << "]";//������
						else str << "],"; //�̾��� ��
						isAvailable = true;
					}
					if (j == polygons.size() - 1) str << "]"; //������
					else str << "],"; //�̾��� ��
				}

				if (i == multiPolygons.size() - 1) str << "]}}"; //���� ��
				else str << "]}}," << endl; //��� �̾��� ��
			}

			if (++iter == tessellated.end()) {
				str << endl;
			}
			else {
				if (isAvailable) str << "," << endl;
			}

			myfile.write(str.str().c_str(), str.str().length());
		}

		myfile.write("]}", 2);
		myfile.close();

	}


private:

	bool debug01 = false;

	double gridUnitGlobal = 100.0;

	//��ü ������ bounding box. �⺻���� EPSG:5179 ��ü ���� ����
	uint32_t BASEX = 740000;
	uint32_t MAXX = 1390000;
	uint32_t BASEY = 1450000;
	uint32_t MAXY = 2070000;

	//���������� ������ �����ϴ� double
	const double EPSILON_DBL = 0.0000001;


	struct GRP_TMP {
		uint32_t idx;
		double posx;
		double posy;
		bool isHead;
		bool isUsed;
	};

	struct GRID_STRG {
		uint32_t beginLoc;
		uint32_t endLoc;
		bool isCCW;
		//bool isUsed;
		vector<uint32_t> segment;
	};

	struct IntersectPt {
		uint32_t gridy;
		double itx;
	};

	struct IDX {
		uint32_t first;
		uint32_t count;
	};

	//������ ���� üũ�ϱ�
	template <typename VECTOR2D>
	bool isPolygonCW(const vector<VECTOR2D>& p) {
		double sum = 0;
		uint32_t vecsize = (uint32_t)p.size();
		for (uint32_t i = 0; i < vecsize; i++) {
			uint32_t j = (i + 1) % vecsize;
			sum += ((double)p[j].x - (double)p[i].x) * ((double)p[j].y + (double)p[i].y);
		}
		return sum > 0;
	}


	//���� ó���� ����, ù ���� ���������� �����ϸ�, ������ �� ����
	//�߰��� �ߺ��� ���鵵 ����
	template <typename VECTOR2D>
	bool removeRedundantPoints(vector<vector<VECTOR2D>>& polygon)
	{
		//uint32_t polygonsize = (uint32_t)polygon.size();
		for (int j = (int)polygon.size() - 1; j >= 0; j--)
		{
			//SPDLOG_DEBUG("loop num : {}", j);
			vector<VECTOR2D>& inner = polygon[j];
			if (inner.size() >= 3) {
				for (int i = (int)inner.size() - 1; i >= 1; i--)
				{
					VECTOR2D& p0 = inner[i - 1];
					VECTOR2D& p1 = inner[i];
					if ((abs(p0.x - p1.x) < EPSILON_DBL) && (abs(p0.y - p1.y) < EPSILON_DBL)) {
						inner.erase(inner.begin() + i);
					}
				}

				VECTOR2D& p0 = inner[0];
				VECTOR2D& p1 = inner[inner.size() - 1];
				if ((abs(p0.x - p1.x) < EPSILON_DBL) && (abs(p0.y - p1.y) < EPSILON_DBL)) {
					inner.erase(inner.begin() + (inner.size() - 1));
				}
			}
			//��ǥ�� 2�� ���϶� �ﰢ�� ������ �ȵǸ� ����ó��
			if (inner.size() <= 2) {
				if (j == 0) return false; //�ٵ� outter�̸� ����ó��
				else polygon.erase(polygon.begin() + j); //�����̸� ����
			}

		}
		return true; //����
	}


	//���� �ٱ���, ccw��, ������ cw�� ����
	//������ ���� ��, �ٱ��� ccw�� �Ǿ�� �ϰ� ������ cw�� �Ǿ�� ��
	template <typename VECTOR2D>
	bool reorderPolygon(vector<vector<VECTOR2D>>& polygon)
	{
		//���� ó���� ����, ù ���� ���������� �����ϸ�, ������ �� ����
		//�߰��� �ߺ��� ���鵵 ����
		bool isOK = removeRedundantPoints(polygon);
		if (!isOK) {
			SPDLOG_ERROR("�ߺ� �� ó�� �� ���� �߻�!!!!!!!!");
			return false;
		}

		uint32_t polygonsize = (uint32_t)polygon.size();
		if (polygonsize < 1) {
			SPDLOG_ERROR("������ ������ 1���� ����. ���� �߻�!!!!!!!!");
			return false;
		}

		bool isOuterPolygonCW = isPolygonCW(polygon[0]);
		if (isOuterPolygonCW) std::reverse(polygon[0].begin(), polygon[0].end());

		for (uint32_t i = 1; i < polygonsize; i++)
		{
			bool isOuterPolygonCW = isPolygonCW(polygon[i]);
			if (!isOuterPolygonCW) std::reverse(polygon[i].begin(), polygon[i].end());
		}

		return true;
	}



	//�׸��忡 �°� ������ �ɰ���. ���۰� ������ �߰��� �����Ѵ�.
	template <typename VECTOR2D>
	vector<VECTOR2D> splitLineByTriangleGrid(VECTOR2D start, VECTOR2D end, uint32_t gridsize,
		vector<IntersectPt>& axis_x,
		bool recordAxisIntersection = false,
		bool addDiagonal = false)
	{
		std::vector<VECTOR2D> segments;
		IntersectPt pt;
		//����
		bool deltaXisZero = (end.x - start.x) == 0;
		const double slope = (end.y - start.y) / (end.x - start.x);

		// Calculate the minimum and maximum x and y values
		double min_x = std::min(start.x, end.x);
		double max_x = std::max(start.x, end.x);
		double min_y = std::min(start.y, end.y);
		double max_y = std::max(start.y, end.y);

		if (addDiagonal) {
			double min_xy = std::min(start.y - start.x, end.y - end.x);
			double max_xy = std::max(start.y - start.x, end.y - end.x);

			//�밢���� �����ϴ� �� �߰�
			int gridStart_xy = min_xy < 0 ? ((int)abs(min_xy) / gridsize) * (-1)
				: ((int)min_xy / gridsize) + 1;
			int gridEnd_xy = max_xy < 0 ? ((int)abs(max_xy) / gridsize) * (-1) - 1
				: ((int)max_xy / gridsize);

			for (int i = gridStart_xy; i <= gridEnd_xy; i++)
			{
				double x = deltaXisZero ? start.x : (start.y - (slope * start.x) - i * gridsize) / (1.0 - slope);
				double y = x + i * gridsize;
				if (y >= min_y && y <= max_y) {
					segments.push_back(VECTOR2D{ x, y });
					//SPDLOG_DEBUG("�� : {}. {}", x, y);
				}
			}
		}

		if (!deltaXisZero)
		{
			uint32_t gridStart_x = (uint32_t)min_x / gridsize + 1;
			uint32_t gridEnd_x = (uint32_t)max_x / gridsize;

			//segments.push_back(start);
			for (uint32_t i = gridStart_x; i <= gridEnd_x; i++) {
				double x = i * gridsize;
				double y = start.y + (x - start.x) * slope;
				if (y >= min_y && y <= max_y) {
					segments.push_back(VECTOR2D{ x, y });
				}
			}
			//segments.push_back(end);
		}

		uint32_t gridStart_y = (uint32_t)min_y / gridsize + 1;
		uint32_t gridEnd_y = (uint32_t)max_y / gridsize;

		for (uint32_t i = gridStart_y; i <= gridEnd_y; i++) {
			double y = i * gridsize;
			double x = deltaXisZero ? start.x : start.x + (y - start.y) * (1.0 / slope);
			if (x >= min_x && x <= max_x) {
				segments.push_back(VECTOR2D{ x, y });
				if (recordAxisIntersection) {
					//������ �׸��� ���� üũ������ �߰��Ѵ�.
					pt.gridy = i * gridsize;
					pt.itx = x;
					axis_x.push_back(pt);
				}
			}
		}

		// ���� ���� ���⿡ ���� �ٸ��� ����
		if (start.x <= end.x) {
			if (start.y <= end.y) std::sort(segments.begin(), segments.end(), [](const VECTOR2D& a, const VECTOR2D& b) {
				return a.x < b.x || (a.x == b.x && a.y < b.y);
				});
			else std::sort(segments.begin(), segments.end(), [](const VECTOR2D& a, const VECTOR2D& b) {
				return a.x < b.x || (a.x == b.x && a.y > b.y);
				});
		}
		else {
			if (start.y <= end.y) std::sort(segments.begin(), segments.end(), [](const VECTOR2D& a, const VECTOR2D& b) {
				return a.x > b.x || (a.x == b.x && a.y < b.y);
				});
			else std::sort(segments.begin(), segments.end(), [](const VECTOR2D& a, const VECTOR2D& b) {
				return a.x > b.x || (a.x == b.x && a.y > b.y);
				});
		}
		segments.erase(std::unique(segments.begin(), segments.end()), segments.end());

		return segments;
	}

	//������ �ϳ��� �޾Ƽ�, �ﰢ�� �׸���� �������ش�.
	template <typename VECTOR2D>
	bool addIntersectedPointsToPolygonByTriangleGrid(vector<vector<VECTOR2D>>& polygons,
		vector<IntersectPt>& axis_x, uint32_t gridunit)
	{
		//�޾ƿ� ������ vector��,
		//ù��° �������� �ٱ�,
		//�ι�° �̻��� �����Ѵٸ� ��� ����
		bool isOK = reorderPolygon(polygons);
		if (!isOK) return false;

		vector<vector<VECTOR2D>> newPolygons; //���� ���� ������

		for (vector<VECTOR2D>& polygon : polygons)
		{
			//newPolygon���� ���� �������� ����, �׸��� ���� ������ ���� ���ԵǾ� �ִ�.
			vector<VECTOR2D> newPolygon;
			newPolygon.push_back(polygon[0]);
			for (uint32_t k = 1; k < polygon.size(); k++) {

				//���� �ϳ��� �׸���� �����Ѵ�. �簢���� �ƴ϶� �ﰢ������ ���ľ� �Ѵ�.
				vector<VECTOR2D> segment = splitLineByTriangleGrid(polygon[k - 1], polygon[k], gridunit, axis_x, true, true);
				segment.push_back(polygon[k]);

				newPolygon.insert(newPolygon.end(), segment.begin(), segment.end());
			}
			//������ ���� ���������Ƿ�, �׺κе� �����ؼ� �߰�
			vector<VECTOR2D> segment = splitLineByTriangleGrid(polygon.back(), polygon.front(), gridunit, axis_x, true, true);
			//segment.push_back(polygon.front());
			newPolygon.insert(newPolygon.end(), segment.begin(), segment.end());
			//�� �������� Ż���ϸ�, newPolygon���� ������ �̱� �����￡ �׸��� ���� ������ ���� ���� ����
			newPolygons.push_back(newPolygon);
		}
		//���� newPolygons�� ������ ���´� �Ȱ�����, �׸���� �����ϴ� ������ �߰��� �������� ��� �ִ�.
		polygons = newPolygons;
		return true;
	}


	template <typename VECTOR2D>
	void getVertexIndirect(const vector<vector<VECTOR2D>>& polygons,
		vector<VECTOR2D>& vertex, vector<IDX>& indirect)
	{
		uint32_t vertexSize = (uint32_t)vertex.size();
		uint32_t countSum = vertexSize; //�̾�޴´�.
		for (const vector<VECTOR2D>& polygon : polygons)
		{
			IDX idx = { countSum, (uint32_t)polygon.size() };
			indirect.push_back(idx);
			countSum += idx.count;
		}
		//ũ�� �����ϰ�

		vertex.resize(countSum);
		uint32_t idx = vertexSize;
		for (const vector<VECTOR2D>& polygon : polygons)
		{
			for (const VECTOR2D& point : polygon)
			{
				vertex[idx] = point;
				idx++;
			}
		}
		//���ؽ��� �ε��� ���� �Ϸ�
	}


	//0 : x����,  1 : y������,  2: �밢��,  3 : ��Ÿ
	template <typename VECTOR2D>
	uint32_t checkPointLoc(VECTOR2D p, uint32_t gridunit)
	{
		uint32_t py = (uint32_t)round(p.y / (double)gridunit) * gridunit;
		if (abs(p.y - py) < EPSILON_DBL) return 0;

		uint32_t px = (uint32_t)round(p.x / (double)gridunit) * gridunit;
		if (abs(p.x - px) < EPSILON_DBL) return 1;

		uint32_t pxy = (uint32_t)round(abs(p.x - p.y) / (double)gridunit) * gridunit;
		if (abs(abs(p.x - p.y) - pxy) < EPSILON_DBL) return 2;

		return 3;
	}


	//�׸���� �߷��� ������ �����̳� �׸���� ��ġ�� ���� ���� �޾Ƶ���
	template <typename VECTOR2D>
	void distributeIdxToGrid(unordered_map<uint32_t, vector<GRID_STRG>>& gridStorage,
		vector<VECTOR2D>& vertex, vector<uint32_t>& idx,
		uint32_t beginLoc, uint32_t endLoc, bool isOuter, uint32_t gridunit, uint32_t xrange)
	{

		GRID_STRG gridstrg = {
			.beginLoc = beginLoc,
			.endLoc = endLoc,
			.isCCW = isOuter,
			.segment = idx
		};

		VECTOR2D pIndicator;
		if (idx.size() > 2) { //���� �� �� �̻��̸� ������ �߰� ��
			pIndicator = vertex[idx[1]];
		}
		else if (idx.size() == 2) { //���� �ΰ����̸� ���� ���
			VECTOR2D p0 = vertex[idx[0]];
			VECTOR2D p1 = vertex[idx[1]];
			double x = (p0.x + p1.x) / 2.0;
			double y = (p0.y + p1.y) / 2.0;
			pIndicator = { x,y };
		}
		else {
			SPDLOG_ERROR("���ҵ� ���п� ���� �ϳ��ۿ� ����!!!!!!!!!");
			exit(1);
		}

		uint32_t gridx = (uint32_t)((pIndicator.x - BASEX) / gridunit);
		uint32_t gridy = (uint32_t)((pIndicator.y - BASEY) / gridunit);
		//�簢�� �׸��带 y = x �������� �������� �� �Ʒ��� 0, ���� 1
		uint32_t xx = gridx * gridunit + BASEX;
		uint32_t yy = gridy * gridunit + BASEY;
		double deltaX = pIndicator.x - xx;
		double deltaY = pIndicator.y - yy;
		uint32_t subidx = deltaX > deltaY ? 0 : 1;
		//SPDLOG_DEBUG("xx : {} , yy : {}, deltaX :{}, deltaY : {}, subidx : {}", xx, yy, deltaX, deltaY, subidx);
		uint32_t gridIndex = (gridy * xrange + gridx) * 2 + subidx;

		gridStorage[gridIndex].push_back(gridstrg);
	}



	template <typename VECTOR2D>
	void distributeSegmentsToGrid(vector<vector<VECTOR2D>>& polygons,
		unordered_map<uint32_t, vector<GRID_STRG>>& gridStorage,
		vector<VECTOR2D>& vertex, unordered_map<uint32_t, uint32_t>& gridPointMap,
		uint32_t gridunit, uint32_t xrange,
		vector<IntersectPt>& axis_x)
	{

		//SPDLOG_DEBUG("indirect �غ�.....");
		vector<IDX> indirect;
		getVertexIndirect(polygons, vertex, indirect);
		//SPDLOG_DEBUG("indirect �Ϸ�..... bbox �غ�");
		//���� ������ �׸��� ������ vertex�� �߰��Ѵ�.

		//�׸��� �ȿ� ���� �߰��ϱ� ���ؼ�, �ϴ� �����Ѵ�.
		std::sort(axis_x.begin(), axis_x.end(), [](const IntersectPt& a, const IntersectPt& b) {
			//if (a.gridy < b.gridy) {
			//	return true;
			//}
			//else if (a.gridy > b.gridy) {
			//	return false;
			//}
			//else {
			//	if (a.itx < b.itx) return true;
			//	else false;
			//}
			return (a.gridy < b.gridy) || ((a.gridy == b.gridy) && (a.itx < b.itx));
			});

		//if (debug01) {
		//	for (IntersectPt& aa : axis_x) SPDLOG_DEBUG("gridy : {}, intersectX : {}", aa.gridy, aa.itx);
		//}


		//�������̶��, �迭 ũ��� ¦���� �Ǿ�� �Ѵ�.
		if (axis_x.size() % 2 == 1) SPDLOG_ERROR("axis_x ¦���� �ƴ�!!!!!!!!!");

		for (uint32_t i = 0; i < axis_x.size(); i += 2)
		{
			IntersectPt& x0 = axis_x[i];
			IntersectPt& x1 = axis_x[i + 1];

			uint32_t xFrom = (uint32_t)(x0.itx / gridunit) * gridunit + gridunit;
			uint32_t xTo = (uint32_t)(x1.itx / gridunit) * gridunit;
			uint32_t y = x0.gridy;
			//SPDLOG_DEBUG("y : {}, xFrom : {} , xTo : {}", y, xFrom, xTo);
			for (uint32_t x = xFrom; x <= xTo; x += gridunit) {
				VECTOR2D xy = { (double)x, (double)y };
				vertex.push_back(xy);
				uint32_t idx_x = (x - BASEX) / gridunit;
				uint32_t idx_y = (y - BASEY) / gridunit;
				uint32_t idx = idx_y * xrange + idx_x; //encode
				//SPDLOG_DEBUG("�߰��� : x: {}, y: {}, idx: {}", x, y, idx);
				gridPointMap[idx] = (uint32_t)vertex.size() - 1; //���߿� ��Ī�� �� �ְ� �ε����� �ִ´�.
			}
		}

		//SPDLOG_DEBUG("�����ϱ�");
		//���� ���鼭 �����Ѵ�.
		bool isOuter = true;
		for (IDX& idx : indirect)
		{
			uint32_t first = idx.first;
			uint32_t last = idx.first + idx.count;

			bool onWriting = false;
			uint32_t beginLoc = 4;
			uint32_t endLoc = 4;
			vector<uint32_t> onWritingIdx;
			vector<uint32_t> skipIdx; //�ʹݿ� �ǳʶٴ� �ε�����

			for (uint32_t i = first; i < last; i++)
			{
				VECTOR2D& p = vertex[i];
				//0 : x����,  1 : y������,  2: �밢��,  3 : ��Ÿ
				uint32_t locCode = checkPointLoc(p, gridunit);

				if (locCode < 3) { //����� �� Ÿ�̹��̸�
					if (onWriting) { //�̹� ������̾�����

						onWritingIdx.push_back(i);
						endLoc = locCode;
						//���⼭ ��� ó��
						distributeIdxToGrid(gridStorage, vertex, onWritingIdx, beginLoc, endLoc, isOuter, gridunit, xrange);

						//���� Ÿ�ֿ̹��� �ε����� �������� ����.
						onWritingIdx.clear();
						onWriting = true;
						onWritingIdx.push_back(i);
						beginLoc = locCode;
					}
					else {
						//������� �ƴϸ� ���� ����. ���� ������ ���� ��ŵ�ϴٰ� ��¥ ó�� �ۿ� ����
						//���� skip �� ���������� ���κ� �߰� �ʿ�
						skipIdx.push_back(i);
						onWriting = true;

						onWritingIdx.push_back(i);
						beginLoc = locCode;
					}
				}
				else { //��� ������̾�� �ϸ�
					if (onWriting) onWritingIdx.push_back(i);
					else skipIdx.push_back(i); //�ʹݿ� ��� �ȵǸ� ���� ������ �ٽ� ������ �ϹǷ� skip�� ����Ѵ�.
				}
			} //for (uint32_t i = first; i < last; i++)

			if (beginLoc == 4)
			{ //beginLoc �� �ʱⰪ�� 4�� ��ϵǾ��ٴ°�, �� ��ü�� �ѹ��� ��� �ȵǾ��ٴ°� �ǹ��Ѵ�.
				//���⼭ ��� ó��. �̹� ��� ��������Ƿ� skipIdx�� ������ �ȴ�.
				distributeIdxToGrid(gridStorage, vertex, skipIdx, 3, 3, isOuter, gridunit, xrange);
			}
			else { //�ܿ� �κ� ��� ó��
				//��ó���� ���� �������� ������ ��ġ�� �ʵ��� �� �ξ���.
				onWritingIdx.insert(onWritingIdx.end(), skipIdx.begin(), skipIdx.end());

				VECTOR2D& p = vertex[skipIdx.back()];
				uint32_t endLoc = checkPointLoc(p, gridunit);
				if (onWritingIdx.size() > 0) distributeIdxToGrid(gridStorage, vertex, onWritingIdx, beginLoc, endLoc, isOuter, gridunit, xrange);
			}


			//ù��° ���� ���ķδ� ��� inner�� ��
			isOuter = false;
		} //for (IDX& idx : indirect)

	}



	bool getGridVertexIndex(uint32_t gridx, uint32_t gridy, uint32_t gridunit, uint32_t xrange,
		unordered_map<uint32_t, uint32_t>& gridPointMap, uint32_t& returnIdx)
	{
		uint32_t idx_x = (gridx - BASEX) / gridunit;
		uint32_t idx_y = (gridy - BASEY) / gridunit;
		uint32_t idx = idx_y * xrange + idx_x; //encode
		if (gridPointMap.contains(idx)) {
			returnIdx = gridPointMap[idx];
			return true;
		}
		else {
			return false;
		}
	}

	bool drawOuterFullTriangle(vector<vector<uint32_t>>& polygon,
		unordered_map<uint32_t, uint32_t>& gridPointMap,
		uint32_t gridx, uint32_t gridy, uint32_t xrange, uint32_t subIdx, uint32_t gridunit)
	{
		vector<uint32_t> triangle;
		uint32_t vertexIndex = 0;
		const uint32_t cornerArr_x = 0b100011; //�ڳ� �迭
		const uint32_t cornerArr_y = 0b101010;
		uint32_t tIdx0 = subIdx * 3;
		for (uint32_t tIdx = tIdx0; tIdx < tIdx0 + 3; tIdx++)
		{
			uint32_t dx = ((cornerArr_x >> tIdx) & 0b01) * gridunit;
			uint32_t dy = ((cornerArr_y >> tIdx) & 0b01) * gridunit;
			bool hasIdx = getGridVertexIndex(
				gridx + dx, gridy + dy,
				gridunit, xrange, gridPointMap, vertexIndex);
			if (hasIdx) {
				triangle.push_back(vertexIndex);
			}
			else {
				SPDLOG_ERROR("ã�� �ִ� �׸��� ��0�� ����!!!!");
				return false;
			}
		}
		triangle.push_back(triangle.front()); //���� ù��° ���� �ٽ� �߰��ؼ� ������
		polygon.push_back(triangle);
		return true;
	}



	template <typename VECTOR2D>
	bool pointInPolygonIdx(VECTOR2D& point, vector<vector<uint32_t>>& polygon, vector<VECTOR2D>& vertex) {

		// ray-casting algorithm based on
		// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
		// http://bl.ocks.org/bycoffe/5575904
		bool intersect;
		bool inside = false;
		for (uint32_t p = 0; p < polygon.size(); p++) {
			uint32_t j = (uint32_t)polygon[p].size() - 1;
			for (uint32_t i = 0; i < polygon[p].size(); j = i++) {
				intersect = ((vertex[polygon[p][i]].y > point.y) != (vertex[polygon[p][j]].y > point.y))
					&& (point.x < (vertex[polygon[p][j]].x - vertex[polygon[p][i]].x)* (point.y - vertex[polygon[p][i]].y) / (vertex[polygon[p][j]].y - vertex[polygon[p][i]].y) + vertex[polygon[p][i]].x);
				if (intersect) inside = !inside;
			}
		}
		return inside;
	}


	template <typename VECTOR2D>
	void constructPolygonFromSegments(unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& tessellated,
		unordered_map<uint32_t, vector<GRID_STRG>>& gridStorage,
		vector<VECTOR2D>& vertex, unordered_map<uint32_t, uint32_t>& gridPointMap,
		uint32_t gridunit, uint32_t xrange)
	{
		for (auto& kv : gridStorage)
		{
			//�׸��� �ϳ��� �ذ��Ѵ�.

			uint32_t grididx = kv.first;
			vector<GRID_STRG>& gridstrg = kv.second;

			//�ε��� �����Ѵ�.
			uint32_t gridxy = grididx / 2;
			uint32_t subIdx = grididx % 2; //�Ʒ� �ﰢ���� true, �� �ﰢ���� false
			uint32_t gridx = (gridxy % xrange) * gridunit + BASEX;
			uint32_t gridy = (gridxy / xrange) * gridunit + BASEY;
			//if (debug01) SPDLOG_DEBUG("================================================================");
			//if (debug01) SPDLOG_DEBUG("gridstorage : {}, {}", gridx, gridy);
			//���庯�� �����ϰ� �ε��� ���
			vector<vector<GRP_TMP>> grp(4); // 0 ����, 1 ������, 2 �밢��, 3 ����


			//�ϴ�, 0, 1, 2, 3 �׷캰�� ����Ѵ�.
			for (uint32_t i = 0; i < gridstrg.size(); i++)
			{
				GRID_STRG& seg = gridstrg[i];
				GRP_TMP grptmp0 = {
					.idx = i,
					.posx = vertex[seg.segment.front()].x, //������
					.posy = vertex[seg.segment.front()].y, //������
					.isHead = true,
					.isUsed = false
				};
				GRP_TMP grptmp1 = {
					.idx = i,
					.posx = vertex[seg.segment.back()].x, //������
					.posy = vertex[seg.segment.back()].y, //������
					.isHead = false,
					.isUsed = false
				};
				grp[seg.beginLoc].push_back(grptmp0);
				if (seg.endLoc != 3) grp[seg.endLoc].push_back(grptmp1); //3���̸� ���� ���̹Ƿ� �ϳ��� �ִ´�.
			}



			//���� �� �ְ� ��Ʈ�Ѵ�.
			if (subIdx == 0) { //�Ʒ��� �ﰢ���̸�
				std::sort(grp[0].begin(), grp[0].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posx < b.posx;
					});
				std::sort(grp[1].begin(), grp[1].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posy < b.posy;
					});
				std::sort(grp[2].begin(), grp[2].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posx > b.posx;
					});

			}
			else { //���� �ﰢ���̸� ��Ʈ ���� �ݴ�
				std::sort(grp[0].begin(), grp[0].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posx > b.posx;
					});
				std::sort(grp[1].begin(), grp[1].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posy > b.posy;
					});
				std::sort(grp[2].begin(), grp[2].end(), [](const GRP_TMP& a, const GRP_TMP& b) {
					return a.posx < b.posx;
					});
			}

			//�߰��ϰ� ��Ʈ�� �� Ȯ�ο�
			if (debug01) {
				uint32_t gi = 0;
				for (vector<GRP_TMP>& g : grp) {
					for (GRP_TMP& t : g) {
						SPDLOG_DEBUG("--------groupIdx : {}, groupSize: {}, segIdx : {}, segIdHead :{}", gi, g.size(), t.idx, t.isHead);
					}
					gi++;
				}
			}

			//�ڿ��� ��ٷ� ����� �ε����� �����.				
			unordered_map<uint32_t, uint32_t> segIdxToVecIdx;
			for (uint32_t i = 0; i <= 2; i++) {
				vector<GRP_TMP>& g = grp[i];
				for (uint32_t j = 0; j < g.size(); j++) {
					uint32_t key = g[j].idx * (3 * 2) + i * 2 + (g[j].isHead ? 0 : 1);
					segIdxToVecIdx[key] = j; //�迭��ȣ�� �ִ´�.
				}
			}
			uint32_t segmentSum = (uint32_t)segIdxToVecIdx.size() / 2; //���߿� while Ż�� üũ��

			//������ �߰��� �� ������ �迭
			//uint32_t cornerArr_x[6] = { gridunit, gridunit, 0,  0, 0, gridunit};
			//uint32_t cornerArr_y[6] = { 0, gridunit, 0,  gridunit, 0, gridunit};
			const uint32_t cornerArr_x = 0b100011; //�ڳ� �迭
			const uint32_t cornerArr_y = 0b101010;


			vector<vector<vector<uint32_t>>> polygons;
			vector<uint32_t> polygonOuter;

			bool isFirst = true;
			uint32_t curGrpIdx = 0;
			uint32_t curIdxInGrp = 0;
			uint32_t usedSegment = 0;
			if (grp[0].size() > 0 || grp[1].size() > 0 || grp[2].size() > 0) {
				while (true)
				{
					//���� üũ&��� ���
					if (debug01) SPDLOG_DEBUG("==curGrpIdx: {}, curIdxInGrp:{}, subIdx:{}", curGrpIdx, curIdxInGrp, subIdx);
					if (isFirst) {
						isFirst = false;
					}
					else {
						//if (curGrpIdx == 0 && curIdxInGrp == 0) break;
					}

					vector<GRP_TMP>& currentGrp = grp[curGrpIdx];
					if (currentGrp.size() > 0) {

						GRP_TMP& g = currentGrp[curIdxInGrp];

						if (g.isHead) {//�Ӹ���
							//SPDLOG_DEBUG("curGrpIdx: {}, curIdxInGrp:{}, isHead! : {}, {}", curGrpIdx, curIdxInGrp, vertex[gridstrg[g.idx].segment[0]].x, vertex[gridstrg[g.idx].segment[0]].y);
							if (!g.isUsed) {
								GRID_STRG& seg = gridstrg[g.idx];
								polygonOuter.insert(polygonOuter.end(), seg.segment.begin(), seg.segment.end());
								usedSegment++;
								//seg.isUsed = true;
								g.isUsed = true;
								curGrpIdx = seg.endLoc;
								uint32_t key = g.idx * (3 * 2) + curGrpIdx * 2 + 1; //������ ã�� ����
								//�ִ��� ������ üũ�ؾ� ��
								if (segIdxToVecIdx.contains(key))
								{
									curIdxInGrp = segIdxToVecIdx[key];
									//SPDLOG_DEBUG("curGrpIdx: {}, curIdxInGrp:{}, �� �߰� �� ����", curGrpIdx, curIdxInGrp);
									//segIdxToVecIdx.erase(key); //�ѹ� ã���� ����
									//continue; //���� ������ ����.
								}
								else { //������ �ȵȴ�.
									SPDLOG_ERROR("segIdxToVecIdx �߸� ����!!!");
									exit(1);
								}
							}
							else { //�Ӹ��ε�, �߰��� ���� ������ �ѹ��� �� ����
								if (polygonOuter.size() > 0) { //���� �߰� ���¸� ��������
									polygonOuter.push_back(polygonOuter.front()); //ù �� �߰��ϰ�
									if (polygonOuter.size() >= 4) { //�ﰢ���� �� ���� 4����. 4�� �̻��̾�� �߰�
										vector<vector<uint32_t>> polygon;
										polygon.push_back(polygonOuter);
										polygons.push_back(polygon); //������ inner�� �����Ƿ� outer�� �߰��ؼ� �����￡ ��Ͻ�Ų��.
									}
									polygonOuter.clear();
									//SPDLOG_DEBUG("segmentSum :{}, usedSegment: {}", segmentSum, usedSegment);
									if (segmentSum == usedSegment) break; //������ �߰��� �����ε�, segment �������� ������ Ż��
								}
								else {//���� �߰� ���°� �ƴϸ�, ��� ���ƾ� ��.
									//�ƹ��͵� ����
								}
							} //if (!g.isUsed) { //�Ӹ���
						} //if (g.isHead)
					}

					//������Ű�鼭 ����� ��
					curIdxInGrp++;
					//SPDLOG_DEBUG("check : curGrp: {},  curIdxInGrp: {}, currentGrp.size():{}", curGrpIdx, curIdxInGrp, grp[curGrpIdx].size());
					if (curIdxInGrp >= grp[curGrpIdx].size())
					{ //���� �׷��� ������ ������,
						if (polygonOuter.size() > 0) { //���� �߰� ���¸�, ������ �߰��ؾ� ��
							uint32_t cornerIdx = subIdx * 3 + curGrpIdx;
							uint32_t dx = ((cornerArr_x >> cornerIdx) & 0b01) * gridunit;
							uint32_t dy = ((cornerArr_y >> cornerIdx) & 0b01) * gridunit;
							uint32_t vertexIndex = 0;
							bool hasIdx = getGridVertexIndex(gridx + dx, gridy + dy, gridunit, xrange, gridPointMap, vertexIndex);

							//SPDLOG_DEBUG("__curGrpIdx: {}, curIdxInGrp:{}", curGrpIdx, curIdxInGrp);
							//SPDLOG_DEBUG("Ȯ���� : x: {}, y: {}, idx: {}", xx, yy, idx);
							if (hasIdx) {
								polygonOuter.push_back(vertexIndex);
							}
							else {
								SPDLOG_ERROR("������ �ȿ� �׸��� ������ �߸� �߰���!!! Ȥ�� ���� ����!!!!!!!!!");
								polygonOuter.clear();
								if (segmentSum == usedSegment) break;
								//exit(1);
							}

						}
						curGrpIdx = (curGrpIdx + 1) % 3; //0 - 1 - 2 ��ȸ
						curIdxInGrp = 0;

					} //if (curIdxInGrp >= currentGrp.size()) 

				} //while (true)
			}

			//������� ������, segment���� �̾����� �������� �������
			//���� inner �� ������ �߰��� ����
			//�ϴ� ccw�� ������ �߰��Ѵ�.
			vector<uint32_t> cws;
			for (uint32_t innerIdx = 0; innerIdx < grp[3].size(); innerIdx++)
			{
				GRP_TMP& inner = grp[3][innerIdx];

				if (gridstrg[inner.idx].isCCW) {

					if (gridstrg[inner.idx].segment.size() >= 3) {
						vector<vector<uint32_t>> polygon;
						polygon.push_back(gridstrg[inner.idx].segment);
						polygon[0].push_back(gridstrg[inner.idx].segment[0]); //ù�� �ٽ� �߰��ؼ� ������ �����
						polygons.push_back(polygon); //������ inner�� �����Ƿ� outer�� �߰��ؼ� �����￡ ��Ͻ�Ų��.
					}
				}
				else { //cw �� ���, 
					cws.push_back(innerIdx); //�ڿ��� ���� ���� ���� ���⼭ �߰�
				}
			}

			//inner �������� �ִµ�, outer�� �ϳ��� ������, ��ü �ﰢ���� �߰��Ѵ�.
			if (cws.size() > 0 && polygons.size() == 0) {
				vector<vector<uint32_t>> polygon;
				bool isOK = drawOuterFullTriangle(polygon, gridPointMap, gridx, gridy, xrange, subIdx, gridunit);
				if (!isOK) continue; //�̹� ���� �ǳʶڴ�.
				polygons.push_back(polygon);
			}

			//�̸� �߰��س��� �ε����� ccw ��������� ��ȯ�Ѵ�.
			for (uint32_t preIdx = 0; preIdx < cws.size(); preIdx++)
			{
				GRP_TMP& inner = grp[3][cws[preIdx]];
				//ccw�� ���� ���� �߰������ϱ� �ǳʶڴ�.
				//if (gridstrg[inner.idx].isCCW) continue;
				//inner �ϳ��� ���鼭 � ������ �������� �Ǻ��Ѵ�.
				for (vector<vector<uint32_t>>& polygon : polygons)
				{
					if (gridstrg[inner.idx].segment.size() >= 3) { //������ 3�� �̻��� ���� �߰�
						VECTOR2D& ptOnSegment = vertex[gridstrg[inner.idx].segment[0]];
						bool isInside = pointInPolygonIdx(ptOnSegment, polygon, vertex);
						if (isInside) {
							polygon.push_back(gridstrg[inner.idx].segment); //���ο� ������ hole�� �߰�
							polygon.back().push_back(gridstrg[inner.idx].segment[0]); //ù�� �ٽ� �߰��ؼ� ������ �����
						}
					}
				}
			}//for (GRP_TMP& inner : grp[3])

			//���� �߰��� �������� ���� �迭�� �ִ´�.
			//���� ������ �̸� �Է��Ѵ�.
			if (!tessellated.contains(gridxy)) {
				vector<vector<vector<uint32_t>>> polygons;
				tessellated[gridxy] = polygons;
			}
			vector<vector<vector<uint32_t>>>& pre = tessellated[gridxy];
			pre.insert(pre.end(), polygons.begin(), polygons.end());

		}// for (auto& kv : gridStorage)

	} //void constructPolygonFromSegments



	template <typename VECTOR2D>
	void fillTrianglesInside(unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& tessellated,
		unordered_map<uint32_t, vector<GRID_STRG>>& gridStorage, //�� ������ ä���� ��� üũ��
		vector<IntersectPt>& axis_x,
		vector<VECTOR2D>& vertex, unordered_map<uint32_t, uint32_t>& gridPointMap,
		uint32_t gridunit, uint32_t xrange)
	{
		if (axis_x.size() < 2) return;
		//������ �� ������ �� �ʿ䰡 �����Ƿ� i < axis_x.size()-2 ������ ����.
		for (uint32_t i = 0; i < axis_x.size() - 2; i += 2)
		{
			IntersectPt& x0 = axis_x[i];
			IntersectPt& x1 = axis_x[i + 1];

			uint32_t xFrom = (uint32_t)(x0.itx / gridunit) * gridunit + gridunit;
			uint32_t xTo = (uint32_t)(x1.itx / gridunit) * gridunit;
			uint32_t gridy = x0.gridy;
			//SPDLOG_DEBUG("y : {}, xFrom : {} , xTo : {}", y, xFrom, xTo);
			for (uint32_t gridx = xFrom; gridx <= xTo; gridx += gridunit) {

				uint32_t gridxx = (uint32_t)((gridx - BASEX) / gridunit);
				uint32_t gridyy = (uint32_t)((gridy - BASEY) / gridunit);
				uint32_t gridxy = (gridyy * xrange + gridxx);
				//SPDLOG_DEBUG("fillTri ���� ��ǥ : {}, {}", gridx, gridy);
				uint32_t gridIndex0 = gridxy * 2 + 0;
				if (!gridStorage.contains(gridIndex0)) {
					//�߰� �� �Ǿ�����, �ﰢ�� �����ؼ� �߰�
					uint32_t subIdx = 0;
					vector<vector<uint32_t>> polygon;
					bool isOK = drawOuterFullTriangle(polygon, gridPointMap, gridx, gridy, xrange, subIdx, gridunit);
					if (isOK) tessellated[gridxy].push_back(polygon);
				}

				uint32_t gridIndex1 = gridxy * 2 + 1;
				if (!gridStorage.contains(gridIndex1)) {
					uint32_t subIdx = 1;
					vector<vector<uint32_t>> polygon;
					bool isOK = drawOuterFullTriangle(polygon, gridPointMap, gridx, gridy, xrange, subIdx, gridunit);
					if (isOK) tessellated[gridxy].push_back(polygon);
				}

			}
		}
	}




}; //class TessPoly

}//namespace vv


#endif 