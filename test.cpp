

#include "tessellatePolygon.hpp"


struct DVEC2 {
	double x;
	double y;

	//std::unique 를 사용하므로 == operator는 반드시 있어야 함
	bool operator==(const DVEC2& other) const {
		return x == other.x && y == other.y;
	}
};



int main() {


	//아웃풋 받을 변수
	unordered_map<uint32_t, vector<vector<vector<uint32_t>>>> tessellatedIndex;
	vector<DVEC2> vertex;

	//그리드사이즈
	double gridsize = 100.0;

	//input으로 넣을 polygon의 bounding box. 충분히 더 넓어도 상관 없다.(계산이 느려지지 않는다)
	uint32_t basex = 740000;
	uint32_t maxx = 1390000;
	uint32_t basey = 1450000;
	uint32_t maxy = 2070000;

	vv::TessPolyGRAT tp;

	//EPSG:5179의 국토 범위가 아닐 경우, init 매개변수에 bounding box를 계산해서 같이 넣어줘야 한다.
	//bounding box의 경우, gridsize로 나누어 떨어지는 숫자가 되어야 한다.
	tp.init(gridsize, basex, maxx, basey, maxy);

	//반드시 구멍을 지닌 폴리곤 형식이어야 한다.
	//구멍이 없다면 바깥 벡터의 원소는 1개가 되어야 한다.
	vector<vector<DVEC2>> polygonWithHole;
	vector<DVEC2> aPolygon = { { 955544, 1953634 }, {955536, 1953593}, {955552, 1953448}, {955565, 1953346}, {955602, 1953180},
		{955821, 1953064}, {955795, 1953058}, {955669, 1952888}, {955671, 1952858}, {955677, 1952830}, {955888, 1952493},
		{955920, 1952453}, {955926, 1952445}, {955989, 1952364}, {955700, 1952294}, {955695, 1952293}, {955662, 1952287},
		{955578, 1952273}, {955562, 1952270}, {955554, 1952268}, {955039, 1952219}, {955014, 1952217}, {954919, 1952210},
		{954865, 1952205}, {954685, 1952217}, {954678, 1952217}, {954662, 1952218}, {954584, 1952237}, {954470, 1952264},
		{954430, 1952272}, {954321, 1952294}, {954317, 1952295}, {953885, 1952337}, {953866, 1952336}, {953818, 1952332},
		{953809, 1953048}, {954048, 1953083}, {954069, 1953079}, {954361, 1953031}, {954467, 1953090}, {954521, 1953120},
		{954573, 1953157}, {954812, 1953233}, {954876, 1953233}, {954906, 1953234}, {954950, 1953774}, {954950, 1953782},
		{954917, 1953935}, {954872, 1954134}, {954799, 1954140}, {954824, 1954353}, {954835, 1954348}, {954912, 1954332},
		{955088, 1954307}, {955176, 1954321}, {955264, 1954291}, {955490, 1953910}, {955503, 1953886}, {955524, 1953821},
		{955548, 1953696}, {955544, 1953634} };
	vector<DVEC2> hole0 = { { 955110, 1953057 }, { 954815, 1952930 }, { 954428, 1952886 }, { 954511, 1952619 },
		{ 955131, 1952575 }, { 955110, 1953057 } };
	vector<DVEC2> hole1 = { { 955267, 1953442 }, { 955250, 1953423 }, { 955267, 1953413 }, { 955293, 1953430 },
		{ 955267, 1953442 } };
	vector<DVEC2> hole2 = { { 955381, 1953339 }, { 955300, 1953352 }, { 955238, 1953300 }, { 955300, 1953233 },
		{ 955469, 1953216 } , { 955381, 1953339 } };
	vector<DVEC2> hole3 = { { 954962, 1953688 }, { 954953, 1953680 }, { 954959, 1953672 }, { 954971, 1953679 },
		{ 954962, 1953688 } };

	polygonWithHole.push_back(aPolygon);
	polygonWithHole.push_back(hole0);
	polygonWithHole.push_back(hole1);
	polygonWithHole.push_back(hole2);
	polygonWithHole.push_back(hole3);
	//테셀레이션
	bool isOK = tp.tessellatePolygon(polygonWithHole, tessellatedIndex, vertex);
	SPDLOG_INFO("그리드 분할 완료");

	//vertex에는 반복되지 않는 정점들이 입력되어 있다.
	//tessellatedIndex 에는 vertex에 대한 index가 있다. 
	//unordered_map의 value 하나하나는 하나의 그리드 단위 공간에 존재하는 multipolygon feature 들이다.(geojson 기준)

	//결과 geojson으로 저장
	string fileName = "d:/temp/result.geojson";
	tp.writeResult(tessellatedIndex, vertex, fileName, (uint32_t)gridsize);
	SPDLOG_INFO("저장 완료 : {}", fileName);
}

