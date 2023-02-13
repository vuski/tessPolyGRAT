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


	//교차점이 추가된 폴리곤을 테셀레이션 한다.
	template <typename VECTOR2D>
	bool tessellatePolygon(vector<vector<VECTOR2D>>& in_polygons,
		unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& out_tessellated,
		vector<VECTOR2D>& out_vertex)
	{
		uint32_t gridunit = (uint32_t)round(gridUnitGlobal);

		//그리드와 교차하는 점들
		vector<IntersectPt> axis_x;

		//홀을 포함한 싱글 폴리곤에 그리드와의 교차점들을 추가한다.
		bool isOK = addIntersectedPointsToPolygonByTriangleGrid(in_polygons, axis_x, gridunit);
		if (!isOK) return false;

		//교차점들이 추가된 폴리곤을 분할한다. 
		//결과로 폴리곤들의 벡터가 반환되어야 한다. 폴리곤들 집합
		//겹치는 점들이 많으므로 버텍스와 인덱스로 관리해야 한다.

		//먼저 폴리곤을 버텍스와 인다이렉트로 변환
		//key는 그리드 인코딩 번호에 삼각형 자른 것
		unordered_map<uint32_t, vector<GRID_STRG>> gridStorage;
		unordered_map<uint32_t, uint32_t> gridPointMap; //좌표 xy인덱스를 넣으면 vertex vector의 인덱스를 반환한다.

		uint32_t xrange = (MAXX - BASEX) / gridunit;
		uint32_t yrange = (MAXY - BASEY) / gridunit;


		//삼각형 그리드 공간에, 폴리곤을 세그먼트로 분할해서 넣는다.
		distributeSegmentsToGrid(in_polygons, gridStorage, out_vertex, gridPointMap, gridunit, xrange, axis_x);
		//SPDLOG_DEBUG("세그먼트 배분 완료");


		//여기까지 오면, 그리드별 저장공간에 폐곡선들이 들어 있음
		//이제 분할된 그리드에 저장된 segment 들을 이어준다.
		//SPDLOG_DEBUG("저장공간에 분할 할당 끝. 이제 구축");
		constructPolygonFromSegments(out_tessellated, gridStorage, out_vertex, gridPointMap, gridunit, xrange);
		//SPDLOG_DEBUG("폴리곤 재구성 완료");


		//마지막으로 중간에 빈 폴리곤을 채운다.
		fillTrianglesInside(out_tessellated, gridStorage, axis_x, out_vertex, gridPointMap, gridunit, xrange);



		return true;
	}


	template <typename VECTOR2D>
	void writeResult(unordered_map<uint32_t, vector<vector<vector<uint32_t>>>>& tessellated,
		vector<VECTOR2D>& vertex, string fileName, uint32_t gridsize)
	{
		uint32_t xrange = (MAXX - BASEX) / gridsize;

		auto startTime = std::chrono::high_resolution_clock::now();
		SPDLOG_DEBUG("tessellated geojson 파일 쓰기 시작.....");
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
						if (k == polygon.size() - 1) str << "]";//끊을때
						else str << "],"; //이어질 때
						isAvailable = true;
					}
					if (j == polygons.size() - 1) str << "]"; //끊을때
					else str << "],"; //이어질 때
				}

				if (i == multiPolygons.size() - 1) str << "]}}"; //끊을 때
				else str << "]}}," << endl; //계속 이어질 때
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

	//전체 데이터 bounding box. 기본값은 EPSG:5179 전체 국토 범위
	uint32_t BASEX = 740000;
	uint32_t MAXX = 1390000;
	uint32_t BASEY = 1450000;
	uint32_t MAXY = 2070000;

	//실질적으로 같음을 정의하는 double
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

	//폴리곤 방향 체크하기
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


	//뒤의 처리를 위해, 첫 점과 마지막점이 동일하면, 마지막 점 제거
	//중간의 중복된 점들도 제거
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
			//좌표가 2개 이하라서 삼각형 구성도 안되면 에러처리
			if (inner.size() <= 2) {
				if (j == 0) return false; //근데 outter이면 에러처리
				else polygon.erase(polygon.begin() + j); //안쪽이면 삭제
			}

		}
		return true; //성공
	}


	//가장 바깥은, ccw로, 안쪽은 cw로 변경
	//폴리곤 정의 상, 바깥은 ccw가 되어야 하고 안쪽은 cw가 되어야 함
	template <typename VECTOR2D>
	bool reorderPolygon(vector<vector<VECTOR2D>>& polygon)
	{
		//뒤의 처리를 위해, 첫 점과 마지막점이 동일하면, 마지막 점 제거
		//중간의 중복된 점들도 제거
		bool isOK = removeRedundantPoints(polygon);
		if (!isOK) {
			SPDLOG_ERROR("중복 점 처리 중 에러 발생!!!!!!!!");
			return false;
		}

		uint32_t polygonsize = (uint32_t)polygon.size();
		if (polygonsize < 1) {
			SPDLOG_ERROR("폴리곤 사이즈 1보다 작음. 에러 발생!!!!!!!!");
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



	//그리드에 맞게 선분을 쪼갠다. 시작과 끝점의 추가는 제외한다.
	template <typename VECTOR2D>
	vector<VECTOR2D> splitLineByTriangleGrid(VECTOR2D start, VECTOR2D end, uint32_t gridsize,
		vector<IntersectPt>& axis_x,
		bool recordAxisIntersection = false,
		bool addDiagonal = false)
	{
		std::vector<VECTOR2D> segments;
		IntersectPt pt;
		//기울기
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

			//대각선과 교차하는 선 추가
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
					//SPDLOG_DEBUG("들어감 : {}. {}", x, y);
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
					//별도로 그리드 범위 체크용으로 추가한다.
					pt.gridy = i * gridsize;
					pt.itx = x;
					axis_x.push_back(pt);
				}
			}
		}

		// 원래 선분 방향에 따라 다르게 정렬
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

	//폴리곤 하나를 받아서, 삼각형 그리드로 분할해준다.
	template <typename VECTOR2D>
	bool addIntersectedPointsToPolygonByTriangleGrid(vector<vector<VECTOR2D>>& polygons,
		vector<IntersectPt>& axis_x, uint32_t gridunit)
	{
		//받아온 폴리곤 vector는,
		//첫번째 폴리곤은 바깥,
		//두번째 이상이 존재한다면 모두 안쪽
		bool isOK = reorderPolygon(polygons);
		if (!isOK) return false;

		vector<vector<VECTOR2D>> newPolygons; //새로 담을 폴리곤

		for (vector<VECTOR2D>& polygon : polygons)
		{
			//newPolygon에는 기존 폴리곤을 따라, 그리드 분할 점들이 같이 포함되어 있다.
			vector<VECTOR2D> newPolygon;
			newPolygon.push_back(polygon[0]);
			for (uint32_t k = 1; k < polygon.size(); k++) {

				//선분 하나를 그리드로 분할한다. 사각형이 아니라 삼각형으로 고쳐야 한다.
				vector<VECTOR2D> segment = splitLineByTriangleGrid(polygon[k - 1], polygon[k], gridunit, axis_x, true, true);
				segment.push_back(polygon[k]);

				newPolygon.insert(newPolygon.end(), segment.begin(), segment.end());
			}
			//마지막 점을 삭제했으므로, 그부분도 감안해서 추가
			vector<VECTOR2D> segment = splitLineByTriangleGrid(polygon.back(), polygon.front(), gridunit, axis_x, true, true);
			//segment.push_back(polygon.front());
			newPolygon.insert(newPolygon.end(), segment.begin(), segment.end());
			//이 루프에서 탈출하면, newPolygon에는 기존의 싱글 폴리곤에 그리드 분할 점들이 같이 찍혀 있음
			newPolygons.push_back(newPolygon);
		}
		//이제 newPolygons에 기존과 형태는 똑같지만, 그리드와 교차하는 점들이 추가된 폴리곤이 들어 있다.
		polygons = newPolygons;
		return true;
	}


	template <typename VECTOR2D>
	void getVertexIndirect(const vector<vector<VECTOR2D>>& polygons,
		vector<VECTOR2D>& vertex, vector<IDX>& indirect)
	{
		uint32_t vertexSize = (uint32_t)vertex.size();
		uint32_t countSum = vertexSize; //이어받는다.
		for (const vector<VECTOR2D>& polygon : polygons)
		{
			IDX idx = { countSum, (uint32_t)polygon.size() };
			indirect.push_back(idx);
			countSum += idx.count;
		}
		//크기 예약하고

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
		//버텍스와 인덱스 생성 완료
	}


	//0 : x수평선,  1 : y수직선,  2: 대각선,  3 : 기타
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


	//그리드로 잘려진 분절된 선분이나 그리드와 겹치지 않은 폐곡선을 받아들임
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
		if (idx.size() > 2) { //점이 두 개 이상이면 적당한 중간 점
			pIndicator = vertex[idx[1]];
		}
		else if (idx.size() == 2) { //점이 두개뿐이면 중점 계산
			VECTOR2D p0 = vertex[idx[0]];
			VECTOR2D p1 = vertex[idx[1]];
			double x = (p0.x + p1.x) / 2.0;
			double y = (p0.y + p1.y) / 2.0;
			pIndicator = { x,y };
		}
		else {
			SPDLOG_ERROR("분할된 선분에 점이 하나밖에 없음!!!!!!!!!");
			exit(1);
		}

		uint32_t gridx = (uint32_t)((pIndicator.x - BASEX) / gridunit);
		uint32_t gridy = (uint32_t)((pIndicator.y - BASEY) / gridunit);
		//사각형 그리드를 y = x 직선으로 분할했을 때 아래가 0, 위가 1
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

		//SPDLOG_DEBUG("indirect 준비.....");
		vector<IDX> indirect;
		getVertexIndirect(polygons, vertex, indirect);
		//SPDLOG_DEBUG("indirect 완료..... bbox 준비");
		//도형 내부의 그리드 점들을 vertex에 추가한다.

		//그리드 안에 점을 추가하기 위해서, 일단 정렬한다.
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


		//정상적이라면, 배열 크기는 짝수가 되어야 한다.
		if (axis_x.size() % 2 == 1) SPDLOG_ERROR("axis_x 짝수가 아님!!!!!!!!!");

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
				//SPDLOG_DEBUG("추가점 : x: {}, y: {}, idx: {}", x, y, idx);
				gridPointMap[idx] = (uint32_t)vertex.size() - 1; //나중에 지칭할 수 있게 인덱스를 넣는다.
			}
		}

		//SPDLOG_DEBUG("분할하기");
		//이제 돌면서 분할한다.
		bool isOuter = true;
		for (IDX& idx : indirect)
		{
			uint32_t first = idx.first;
			uint32_t last = idx.first + idx.count;

			bool onWriting = false;
			uint32_t beginLoc = 4;
			uint32_t endLoc = 4;
			vector<uint32_t> onWritingIdx;
			vector<uint32_t> skipIdx; //초반에 건너뛰는 인덱스들

			for (uint32_t i = first; i < last; i++)
			{
				VECTOR2D& p = vertex[i];
				//0 : x수평선,  1 : y수직선,  2: 대각선,  3 : 기타
				uint32_t locCode = checkPointLoc(p, gridunit);

				if (locCode < 3) { //끊어야 할 타이밍이면
					if (onWriting) { //이미 기록중이었으면

						onWritingIdx.push_back(i);
						endLoc = locCode;
						//여기서 기록 처리
						distributeIdxToGrid(gridStorage, vertex, onWritingIdx, beginLoc, endLoc, isOuter, gridunit, xrange);

						//끊는 타이밍에는 인덱스가 이중으로 들어간다.
						onWritingIdx.clear();
						onWriting = true;
						onWritingIdx.push_back(i);
						beginLoc = locCode;
					}
					else {
						//기록중이 아니면 새로 시작. 여기 들어오는 경우는 스킵하다가 진짜 처음 밖에 없음
						//따라서 skip 의 마지막에도 경계부분 추가 필요
						skipIdx.push_back(i);
						onWriting = true;

						onWritingIdx.push_back(i);
						beginLoc = locCode;
					}
				}
				else { //계속 기록중이어야 하면
					if (onWriting) onWritingIdx.push_back(i);
					else skipIdx.push_back(i); //초반에 기록 안되면 루프 끝나고 다시 돌려야 하므로 skip에 기록한다.
				}
			} //for (uint32_t i = first; i < last; i++)

			if (beginLoc == 4)
			{ //beginLoc 이 초기값인 4로 기록되었다는건, 폐곡선 전체가 한번도 기록 안되었다는걸 의미한다.
				//여기서 기록 처리. 이미 모두 들어있으므로 skipIdx를 넣으면 된다.
				distributeIdxToGrid(gridStorage, vertex, skipIdx, 3, 3, isOuter, gridunit, xrange);
			}
			else { //잔여 부분 기록 처리
				//전처리를 통해 시작점과 끝점이 겹치지 않도록 해 두었다.
				onWritingIdx.insert(onWritingIdx.end(), skipIdx.begin(), skipIdx.end());

				VECTOR2D& p = vertex[skipIdx.back()];
				uint32_t endLoc = checkPointLoc(p, gridunit);
				if (onWritingIdx.size() > 0) distributeIdxToGrid(gridStorage, vertex, onWritingIdx, beginLoc, endLoc, isOuter, gridunit, xrange);
			}


			//첫번째 바퀴 이후로는 계속 inner가 됨
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
		const uint32_t cornerArr_x = 0b100011; //코너 배열
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
				SPDLOG_ERROR("찾고 있는 그리드 점0이 없음!!!!");
				return false;
			}
		}
		triangle.push_back(triangle.front()); //가장 첫번째 점을 다시 추가해서 마무리
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
			//그리드 하나씩 해결한다.

			uint32_t grididx = kv.first;
			vector<GRID_STRG>& gridstrg = kv.second;

			//인덱스 복원한다.
			uint32_t gridxy = grididx / 2;
			uint32_t subIdx = grididx % 2; //아래 삼각형은 true, 위 삼각형은 false
			uint32_t gridx = (gridxy % xrange) * gridunit + BASEX;
			uint32_t gridy = (gridxy / xrange) * gridunit + BASEY;
			//if (debug01) SPDLOG_DEBUG("================================================================");
			//if (debug01) SPDLOG_DEBUG("gridstorage : {}, {}", gridx, gridy);
			//저장변수 셋팅하고 인덱스 배분
			vector<vector<GRP_TMP>> grp(4); // 0 수평선, 1 수직선, 2 대각선, 3 독립


			//일단, 0, 1, 2, 3 그룹별로 배분한다.
			for (uint32_t i = 0; i < gridstrg.size(); i++)
			{
				GRID_STRG& seg = gridstrg[i];
				GRP_TMP grptmp0 = {
					.idx = i,
					.posx = vertex[seg.segment.front()].x, //기준점
					.posy = vertex[seg.segment.front()].y, //기준점
					.isHead = true,
					.isUsed = false
				};
				GRP_TMP grptmp1 = {
					.idx = i,
					.posx = vertex[seg.segment.back()].x, //기준점
					.posy = vertex[seg.segment.back()].y, //기준점
					.isHead = false,
					.isUsed = false
				};
				grp[seg.beginLoc].push_back(grptmp0);
				if (seg.endLoc != 3) grp[seg.endLoc].push_back(grptmp1); //3번이면 단일 폐곡선이므로 하나만 넣는다.
			}



			//이을 수 있게 소트한다.
			if (subIdx == 0) { //아래쪽 삼각형이면
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
			else { //위쪽 삼각형이면 소트 순서 반대
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

			//추가하고 소트한 점 확인용
			if (debug01) {
				uint32_t gi = 0;
				for (vector<GRP_TMP>& g : grp) {
					for (GRP_TMP& t : g) {
						SPDLOG_DEBUG("--------groupIdx : {}, groupSize: {}, segIdx : {}, segIdHead :{}", gi, g.size(), t.idx, t.isHead);
					}
					gi++;
				}
			}

			//뒤에서 곧바로 사용할 인덱스를 만든다.				
			unordered_map<uint32_t, uint32_t> segIdxToVecIdx;
			for (uint32_t i = 0; i <= 2; i++) {
				vector<GRP_TMP>& g = grp[i];
				for (uint32_t j = 0; j < g.size(); j++) {
					uint32_t key = g[j].idx * (3 * 2) + i * 2 + (g[j].isHead ? 0 : 1);
					segIdxToVecIdx[key] = j; //배열번호를 넣는다.
				}
			}
			uint32_t segmentSum = (uint32_t)segIdxToVecIdx.size() / 2; //나중에 while 탈출 체크용

			//꼭지점 추가할 때 참조할 배열
			//uint32_t cornerArr_x[6] = { gridunit, gridunit, 0,  0, 0, gridunit};
			//uint32_t cornerArr_y[6] = { 0, gridunit, 0,  gridunit, 0, gridunit};
			const uint32_t cornerArr_x = 0b100011; //코너 배열
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
					//현재 체크&등록 대상
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

						if (g.isHead) {//머리면
							//SPDLOG_DEBUG("curGrpIdx: {}, curIdxInGrp:{}, isHead! : {}, {}", curGrpIdx, curIdxInGrp, vertex[gridstrg[g.idx].segment[0]].x, vertex[gridstrg[g.idx].segment[0]].y);
							if (!g.isUsed) {
								GRID_STRG& seg = gridstrg[g.idx];
								polygonOuter.insert(polygonOuter.end(), seg.segment.begin(), seg.segment.end());
								usedSegment++;
								//seg.isUsed = true;
								g.isUsed = true;
								curGrpIdx = seg.endLoc;
								uint32_t key = g.idx * (3 * 2) + curGrpIdx * 2 + 1; //꼬리를 찾는 것임
								//있는지 없는지 체크해야 함
								if (segIdxToVecIdx.contains(key))
								{
									curIdxInGrp = segIdxToVecIdx[key];
									//SPDLOG_DEBUG("curGrpIdx: {}, curIdxInGrp:{}, 선 추가 후 변경", curGrpIdx, curIdxInGrp);
									//segIdxToVecIdx.erase(key); //한번 찾으면 삭제
									//continue; //다음 루프로 간다.
								}
								else { //없으면 안된다.
									SPDLOG_ERROR("segIdxToVecIdx 잘못 만듬!!!");
									exit(1);
								}
							}
							else { //머리인데, 추가된 적이 있으면 한바퀴 돈 것임
								if (polygonOuter.size() > 0) { //현재 추가 상태면 마지막임
									polygonOuter.push_back(polygonOuter.front()); //첫 점 추가하고
									if (polygonOuter.size() >= 4) { //삼각형일 때 점이 4개임. 4개 이상이어야 추가
										vector<vector<uint32_t>> polygon;
										polygon.push_back(polygonOuter);
										polygons.push_back(polygon); //아직은 inner가 없으므로 outer만 추가해서 폴리곤에 등록시킨다.
									}
									polygonOuter.clear();
									//SPDLOG_DEBUG("segmentSum :{}, usedSegment: {}", segmentSum, usedSegment);
									if (segmentSum == usedSegment) break; //폴리곤 추가한 상태인데, segment 남아있지 않으면 탈출
								}
								else {//현재 추가 상태가 아니면, 계속 돌아야 함.
									//아무것도 안함
								}
							} //if (!g.isUsed) { //머리면
						} //if (g.isHead)
					}

					//증가시키면서 맞춰야 함
					curIdxInGrp++;
					//SPDLOG_DEBUG("check : curGrp: {},  curIdxInGrp: {}, currentGrp.size():{}", curGrpIdx, curIdxInGrp, grp[curGrpIdx].size());
					if (curIdxInGrp >= grp[curGrpIdx].size())
					{ //현재 그룹의 끝까지 갔으면,
						if (polygonOuter.size() > 0) { //현재 추가 상태면, 꼭지점 추가해야 함
							uint32_t cornerIdx = subIdx * 3 + curGrpIdx;
							uint32_t dx = ((cornerArr_x >> cornerIdx) & 0b01) * gridunit;
							uint32_t dy = ((cornerArr_y >> cornerIdx) & 0b01) * gridunit;
							uint32_t vertexIndex = 0;
							bool hasIdx = getGridVertexIndex(gridx + dx, gridy + dy, gridunit, xrange, gridPointMap, vertexIndex);

							//SPDLOG_DEBUG("__curGrpIdx: {}, curIdxInGrp:{}", curGrpIdx, curIdxInGrp);
							//SPDLOG_DEBUG("확인점 : x: {}, y: {}, idx: {}", xx, yy, idx);
							if (hasIdx) {
								polygonOuter.push_back(vertexIndex);
							}
							else {
								SPDLOG_ERROR("폴리곤 안에 그리드 점들을 잘못 추가함!!! 혹은 위상 에러!!!!!!!!!");
								polygonOuter.clear();
								if (segmentSum == usedSegment) break;
								//exit(1);
							}

						}
						curGrpIdx = (curGrpIdx + 1) % 3; //0 - 1 - 2 순회
						curIdxInGrp = 0;

					} //if (curIdxInGrp >= currentGrp.size()) 

				} //while (true)
			}

			//여기까지 왔으면, segment들이 이어져서 폴리곤을 만들었음
			//이제 inner 가 있으면 추가할 차례
			//일단 ccw가 있으면 추가한다.
			vector<uint32_t> cws;
			for (uint32_t innerIdx = 0; innerIdx < grp[3].size(); innerIdx++)
			{
				GRP_TMP& inner = grp[3][innerIdx];

				if (gridstrg[inner.idx].isCCW) {

					if (gridstrg[inner.idx].segment.size() >= 3) {
						vector<vector<uint32_t>> polygon;
						polygon.push_back(gridstrg[inner.idx].segment);
						polygon[0].push_back(gridstrg[inner.idx].segment[0]); //첫점 다시 추가해서 폴리곤 만들기
						polygons.push_back(polygon); //아직은 inner가 없으므로 outer만 추가해서 폴리곤에 등록시킨다.
					}
				}
				else { //cw 일 경우, 
					cws.push_back(innerIdx); //뒤에서 루프 돌기 위해 여기서 추가
				}
			}

			//inner 폴리곤은 있는데, outer가 하나도 없으면, 전체 삼각형을 추가한다.
			if (cws.size() > 0 && polygons.size() == 0) {
				vector<vector<uint32_t>> polygon;
				bool isOK = drawOuterFullTriangle(polygon, gridPointMap, gridx, gridy, xrange, subIdx, gridunit);
				if (!isOK) continue; //이번 판은 건너뛴다.
				polygons.push_back(polygon);
			}

			//미리 추가해놓은 인덱스로 ccw 폴리곤들을 순환한다.
			for (uint32_t preIdx = 0; preIdx < cws.size(); preIdx++)
			{
				GRP_TMP& inner = grp[3][cws[preIdx]];
				//ccw는 조금 전에 추가했으니까 건너뛴다.
				//if (gridstrg[inner.idx].isCCW) continue;
				//inner 하나씩 돌면서 어떤 폴리곤 내부인지 판별한다.
				for (vector<vector<uint32_t>>& polygon : polygons)
				{
					if (gridstrg[inner.idx].segment.size() >= 3) { //정점이 3개 이상일 때만 추가
						VECTOR2D& ptOnSegment = vertex[gridstrg[inner.idx].segment[0]];
						bool isInside = pointInPolygonIdx(ptOnSegment, polygon, vertex);
						if (isInside) {
							polygon.push_back(gridstrg[inner.idx].segment); //내부에 있으면 hole로 추가
							polygon.back().push_back(gridstrg[inner.idx].segment[0]); //첫점 다시 추가해서 폴리곤 만들기
						}
					}
				}
			}//for (GRP_TMP& inner : grp[3])

			//이제 추가한 폴리곤을 기존 배열에 넣는다.
			//아직 없으면 미리 입력한다.
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
		unordered_map<uint32_t, vector<GRID_STRG>>& gridStorage, //이 변수는 채워진 경계 체크용
		vector<IntersectPt>& axis_x,
		vector<VECTOR2D>& vertex, unordered_map<uint32_t, uint32_t>& gridPointMap,
		uint32_t gridunit, uint32_t xrange)
	{
		if (axis_x.size() < 2) return;
		//마지막 맨 윗줄은 할 필요가 없으므로 i < axis_x.size()-2 까지만 돈다.
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
				//SPDLOG_DEBUG("fillTri 점검 좌표 : {}, {}", gridx, gridy);
				uint32_t gridIndex0 = gridxy * 2 + 0;
				if (!gridStorage.contains(gridIndex0)) {
					//추가 안 되었으면, 삼각형 구성해서 추가
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