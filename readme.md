## 폴리곤을 직각삼각형 격자에 맞춰 분할하기
## Tessellating a Polygon into Grids of Right-Angled Triangles

- tessellatePolygon.hpp 는 폴리곤을 직각삼각형 격자로 분할합니다.
- 아래와 같은 폴리곤을 넣으면
![](/image/input1.png)

- 아래와 같은 폴리곤을 만들어줍니다.
![](/image/output1.png)

- 반복 사용되는 점들이 많으므로 vertex와 index 형식으로 결과를 내놓습니다.

- test.cpp 에는 간단한 사용 예시가 있습니다. 

- writeResult 함수는 결과물을 geojson 으로 저장합니다. test.cpp 의 마지막 저장 함수 내용을 보시면 출력 데이터 구조를 이해하실 수 있습니다.

- EPSG:5179 좌표계의 대한민국 영토 내부의 폴리곤이 기본 대상으로 설정되어 있으나, init 함수에서 bounding box값을 다르게 설정하면 다른 영역의 도형도 분할 가능합니다.

- tessellatePolygon 함수는 하나의 폴리곤에 대해서만 작동하며, 다수의 폴리곤으로 이루어진 경우에는 루프를 구성해서 전후 처리를 해주어야 합니다.

- 화면 출력 용도로 spdlog 라이브러리를 사용합니다. 함께 동봉되어 있습니다.
