на текущий момент есть следующие теоремы:
	- nz6, nz8
	- Z_6 connectivity
	- 7c4c
	- 10c6c (и вообще ckc, где k >= 4)
	- o11c6c
	- o14c8c (и вообще ockc, где k >= 12)


на текущий момент есть следующие сильные не следующие друг из друга гипотезы:
	- 2bm с dominating circuit => nz5, 5cdc
	- 3bm => nz5, o7c4c
	- petersen colouring => 5cdc, 6c4c, 10c6c (который оптимален кажется только для графа петерсена; хотя у него и o10c6c есть)
	- 5-face-colourable circular embedding => o5cdc => nz5, 5cdc
	- o6c4c => 6c4c => 22/15
	- nz5 на сфере + unit vector flows => nz5
	- 21/15 => cdc
	- Z_5 connectivity => nz5
	
	- 333-flows
	- oriented 334-flows => o10c4c
	- oriented 244-flows => o10c4c
	- hoffman-ostenhoff
	- o9c6c (кроме графа петерсена)
	
	- faithful circuit cover (тут вроде бы и 2 варианта, но я ж так понимаю, что 1) => cdc
	- small oriented cdc => small cdc [и oppdc ещё] => cdc
	
	- oriented shortest circuit cover? (всё же кажется что это 21/15 опять же и тогда из 21/15 => ocdc в принципе)

16 гипотез
хотя и мб, что 21/15, oriented shortest circuit cover и socdc все про одно и то же
но в любом случае тогда имеем по-крайней мере 14 гипотез

вообще говоря
faithful circuit cover, скажем, испольует 3 понятия
	admissible, eulerian и cut
	и типа
	это уже намёк на какую-то структуру
	которая могла бы пригодиться в других понятиях
в Z_5 connectivity тоже есть своё понятие: boundary
в 2bm - dominating circuit и разбиение многих оставшихся рёбер на 2 множества
в 3bm - набор circuit, которые в итоге dominating и разбиение оставшихся рёбер на 3 множества
в 2bm и 3bm есть ignored вершины, не попавшие в dominating circuit(s)
в hoffman-ostenhoff есть spanning tree
	в nz6/Z_6 connectivity он вроде тоже есть?
	в Z_6 есть TC3, и в этих TC3 вроде как есть 3 spanning tree (double cover)
	в nz6 есть разбиение на 2 потока, nz3 и nz2
в petersen colouring есть вложение в граф petersen'а, есть poor, есть rich рёбра
	в 21/15 кстати интересно, вроде там нет такого сведения, но граф петерсена участвует в доказательстве 21/15 => cdc
	в o6c4c (да и в 6c4c) тоже есть poor и rich рёбра
в o6c4c есть oriented вершины
в o5cdc, o6c4c, oriented 334-flows, oriented 244-flows и small oriented cdc есть ориентация на циклах
	возможно что в 21/15 она тоже есть
в 333-flows, oriented 334-flows, oriented 244-flows хотелось бы какого-то усиления
	в плане что там же 3 подграфа получается с этими потоками
	хочется от них чего-то спросить, какого-то ограничения
	ну то есть как не самый лучший тут пример - oriented 244-flows
	2 - это цикл, чётной длины; половина рёбер накрывается одним из 4-flow графом, вторая половина - другим
в nz5 на сфере + unit vector flows есть двумерная сфера
	а может это и не сфера а проективная плоскость из-за того, что сумма на противоположных точках нулевая
в nz5 на сфере + unit vector flows, 2bm, 3bm, o5cdc, Z_5 connectivity есть nz5
	более того, в 2bm и o5cdc есть 5cdc и nz5
	
- TODO: а что у матроидов?
	- у них точно есть аналог cdc (но не знаю про k-cdc и oriented cdc и small cdc и shortest cdc и strong cdc)
	- c4c? c6c?
	- насколько я понимаю, у матроидов бывают всякие разные k-nowhere zero flows, и есть история/связь с shortest cycle cover
	- dominating circuit? bipartizing matchings?
	- group connectivity?
	- k-flow matroid double cover?
	- faithful circuit cover?
	- hoffman-ostenhoff
	
	- signed матроиды?
		- signed-graphic matroid (frame matroid)
		- extended lift matroid

	- ну заодно можно было бы попробовать поизучать и другие виды матроидов, типа delta-matroids
	
- signed графы:
	должны быть какие-то signed аналоги примерно у всех гипотез
	но непонятно как перебрать все графы и все различные signed-подмножества-рёбер в каждом графе
	(и опять же - какие именно брать графы, потому что там есть куча контрпримеров уже хотя бы для signed cdc)
	signed nz6
	signed [o]6cdc
	signed [o]6c4c
	signed [o]9?c6c (и отдельно signed [o]10?c6c для графа Петерсена)
	signed Z_k connectivity (Z_6?)
	signed dominating circuit
	signed hoffman-ostenhoff? (скажем, чтоб цикл был balanced или может даже просто signed цикл)
		A weaker question is whether
		a 2-connected signed graph with an even number of negative edges has a circuit double cover or not? The
		answer to this question is negative, even for 3-connected cubic signed graph. The signed graph in Figure 1
		has no circuit double cover.
		As many 2-edge-connected signed graphs have no circuit double covers, it is interesting to ask, is there
		an integer k such that every 2-connected flow-admissible signed graph (G, σ) has a circuit k-cover?

		кстати
		The negativeness of a signed graph (G, σ) is the smallest number of negative edges over all equivalent
		signatures of σ, denoted by ǫ(G, σ). M´aˇcajov´a and Skoviera [16] proved that a 2-edge-connected signed ˇ
		graph is flow-admissible if and only if ǫ(G, σ) != 1. Combining it with Boucet’s result [3] that a signed
		graph with a circuit cover if and only if it is flow-admissible, the following observation holds.
		Observation 2.2. Let (G, σ) be a 2-edge-connected signed graph. Then (G, σ) has a circuit cover if and
		only if ǫ(G, σ) != 1.

нужно для всех мелких снарков (ну хотя бы для начала для снарков на 10,18,20 вершинах; итого 1+2+6=9 снарков) описать по возможности все:
	- o5cdc
	- o6c4c
	- 5cdc
	- 6c4c
	- 2bm
	- 3bm
	- dominating circuit
	- petersen colouring
	- hoffman-ostenhoff
	- nz5
	- nz-k polynomial
	- Z_5 connectivity
	- shortest 4-cycle cover aka 21/15 + shortest oriented circuit cover
	- shortest 3-cycle cover aka 22/15
	- 333-flows
	- oriented 334-flows
	- oriented 244-flows
	- 4-perfect matching covering (тут вопрос скорее, есть оно или нет его)
	- unit vector flows
	- o9c6c (причём наверняка ещё найдутся контрпримеры)
	- 9c6c
	- oriented shortest circuit cover
	- dot product decompositions
	- small [oriented] cdc (хотя тут наверняка всегда решения мелкие)
	- faithful circuit cover
	- более маргинальные вещи: [oriented] 2223-flows, antisymmetric flow, 3-edge connectedness, ppdc, oppdc, eppdc, edp- pp- spp- flows, orientation для perfect matchings в o6c4c
и потом поискать точки соприкосновения
и ещё может где-то отдельно поизучать
	- nz3 for any graph degree sequence

скажем,
10.05g1:
	- o5cdc: 96555 - на торе - это ещё и 3bm с 1 ignored
	- 5cdc: добавляется 86655 - на бутылке клейна, то есть неориентируемая поверхность
	- o6c4c - 55;55;55;55;55;55, где 3 oriented вершины (наверняка это тоже куда-то вкладывается на какой-то аналог поверхности, но я не придумал конструкцию) - эти вершины являются соседями одной единственной вершины в графе; все рёбра RICH
	- 6с4с - то же решение
	- 2bm, 3bm, dominating circuit: все строятся из цикла в 9 вершин
	- stronger petersen colouring, petersen colouring: 1 решение, где все рёбра RICH
	- hoffman-ostenhoff: 2 решения - 6 + дерево; 5 + дерево + ребро
	- nz5: их 2400 штук, если не учитывать симметрий (видимо если учитывать, то получится 20 решений)
	- nz-k polynomial: 2400, 19080, 85080(, 278880 и к счастью это число совпало с вычислением по полиному)
		вообще интересно, как это компактно описать
	- shortest 4-cycle cover (aka 21/15): 6555
	- shortest 3-cycle cover (aka 22/15):
		Every bridgeless cubic graph with m edges has a 3-even subgraph cover with total length at most 22/15m
		10+6+6? таких нет решений (типа осталось 5 рёбер; но каждый из циклов длины 6 может накрыть только 2 ребра из 5)
		8+8+6? смотрю по geogebra, 10g1: 08372916, 15429806, 045376
		получилось
		как они пересекаются:
			1-2: 08, 29, 16, 06
			1-3: 06, 37
			2-3: 45, 06
		хм, тут есть ребро, которое накрыто 3 раза
		может я могу найти решение,где нет такого спецэффекта?
		допустим есть циклы 08372916, 045376
		нужен ещё цикл длины 8
		в котором нет рёбер 06, 37, но есть рёбра 15, 24, 89
		не, нет такого цикла
		
		допустим цикл длины 6 расположен по-другому (а это всего 1 неизоморфно различный способ)
		скажем 160835
		но тогда одним циклом не накроешь все 3 ребра у вершины 4
		ну всё, значит от этого спецэффекта не избавиться
		то есть: 3 раза накрыто ребро 06; 2 раза - 08, 29, 16, 37, 45 (кстати, это perfect matching); 1 раз все остальные
		3*1 + 2*5 + 1*9 = 22
	- 233-flows, 234-flows, 235-flows (и вероятно вообще 23k-flows): ничего нет
	- 4-perfect matching covering: для графа Петерсена их не существует
	- unit vector flow:
		ребро графа петерсена - это пара противоположных вершин icosidodecahedron'а
		пара - потому что ребро можно повернуть в 2 разных стороны
	- o9c6c, 9c6c - не существует таких (доказал)
		10c6c существуют, и довольно красивые: скажем, 10 циклов длины 9
		o10c6c тоже существуют, но там уже логики я не нашёл (кроме того, что должны присутствовать всегда циклы длины 5; а вообще там всякие решения бывают, даже что 2 слоя противоположно направлены)
	- dot product decompositions: нет таких, потому что это самый простой снарк
	
	- TODO: Z_5 connectivity
	- TODO: 333-flows
	- TODO: oriented 334-flows
	- TODO: oriented 244-flows

18.05g1:
	- o5cdc: как минимум 24 решения (хотя некоторые друг на друга похожи - кажется, что можно перебрасывать циклы между слоями)
		2 из этих 24 решений содержат dominating circuit
			17;5+6;5;6+6;9; 
			17;5+6;5+6;6;9;
	- 5cdc: ещё как минимум >= 477 решений
	- stronger petersen colouring, petersen colouring: 1 неизоморфное решение, где есть 2 пути по 4 poor ребра

	o5cdc
	глянем решения с циклом длиной 17
	что попадает в ignored: 16, 17
	и всё
	в stronger petersen colouring: эти 2 вершины всегда являются центрами путей из poor рёбер
	в 2/3bm: в ignored попадает всё что угодно (при длине цикла 17)

	
18.05g2:
	- o5cdc: как минимум 54 решения
		2 из этих 24 решений содержат dominating circuit
			17;5+6;5;6+6;9; 
			17;5+6;5+6;6;9;
	- 5cdc: ещё как минимум >= 665 решений
	- stronger petersen colouring, petersen colouring: 1 неизоморфное решение, где есть 2 пути по 4 poor ребра
	- nz5: 226416 решений
	- nz-k polynomial: 226416, 7081284, ...
	
	o5cdc
	ignored для циклов длины 17: 8, 11
	и тоже всё
	в stronger petersen colouring: эти 2 вершины всегда являются концами путей из poor рёбер

20.05g1:
0 : 10 4 14
1 : 9 18 6
2 : 15 4 11
3 : 8 19 7
4 : 0 2 5
5 : 4 12 16
6 : 1 10 7
7 : 3 6 15
8 : 3 14 9
9 : 1 8 11
10 : 0 6 13
11 : 2 9 13
12 : 5 18 13
13 : 10 11 12
14 : 0 8 17
15 : 2 7 17
16 : 5 19 17
17 : 14 15 16
18 : 1 12 19
19 : 3 16 18
	у этого графа особенно выделяется один цикл длины 5 (собственно он в нём один такой)
		'12 18 19 16 5'
	- o5cdc - есть решения с этим циклом
	- 2bm: есть всякие решения, включающие в себя этот цикл
		есть даже такое:
			circuit lengths: 10 5
			circuits:
			10 6 7 15 2 11 9 8 14 0
			12 18 19 16 5
	- stronger petersen colouring: есть решения, где 9 poor рёбер - 2 пути по 2 ребра + этот цикл
		если не учитывать симметрии, то
		2 решения с 5 poor рёбрами
		10 - с 7
		5 - с 9 (во всех есть этот цикл)
	- hoffman-ostenhoff:
		сначала неправильные недорешения:
			есть 2 интересных недорешения, связанные с этим циклом stronger petersen colouring и прочим
			а именно - в качестве цикла возьмём этот самый цикл
			в качестве matching - один из двух наборов, которые бывают poor:
				(0,14) (2,11) (6,10) (7,15) (8,9) 
				(0,10) (2,15) (6,7) (8,14) (9,11) 
			и в качестве почтидерева - всё остальное
		можно построить правильное решение: в качестве цикла взять '0 14 8 9 11 2 15 7 6 10'
		в качестве matching - любое ребро из цикла длины 5
		в качестве дерева - всё остальное
		т. е. получаем 10+1+19
		
	- o5cdc:
		нет нигде цикла длиной 19
		единственный цикл длины 5 во всех решениях - это '12 18 19 16 5'
			это вообще единственный цикл длины 5 в графе
		18:
			8,15; 9,12; 7,14; 0,16; 2,12; 0,12;
		в strong (? TODO) petersen нет следующих вершин (они никогда не граничат с poor рёбрами):
			1, 3, 4, 13, 17
			
		кстати прикольно
		есть один из наборов слабых рёбер - (0,14) (2,11) (6,10) (7,15) (8,9) ; есть цикл '12 18 19 16 5'; и есть оставшиеся вершины 1, 3, 4, 13, 17



20.05g2:
	- petersen colouring: 5 решений (2 с 9 poor, 2 с 11 poor, 1 с 13 poor)
		видимо 3 неизоморфных решения
		2 из них - не stronger (где 11 или 13 poor рёбер)
		рёбра 2-7, 2-9 всегда poor
		рёбра, которые всегда rich: 0-12, 6-7, 9-11
		вершины, которые всегда соседи с poor: 5
		все вершины так или иначе бывают соседями с poor рёбрами

20.05g3:
	- petersen colouring: 10 решений, неизоморфных ~ 2
	всегда poor рёбра: 1-11
	всегда rich рёбра: 1-14, 11-17, 18-19
	вершины, всегда соседние с poor рёбрами: 1, 11
	вершины, которые никогда не соседние с poor рёбрами: -
	все решения - "stronger" petersen colouring:
		- пути: 4, циклы: 7
		- пути: 2+3+4
	
20.05g4:
	- petersen colouring: 7 решений, неизоморфных ~ 4
	всегда poor рёбра: 1-11
	всегда rich рёбра: 10-12
	вершины, всегда соседние с poor рёбрами: 1, 11
	вершины, которые никогда не соседние с poor рёбрами: -
	есть не stronger решения
	"stronger" petersen colouring:
		- пути: 4, циклы: 7
		- пути: 2+3+4
	
20.05g5:
	- petersen colouring: 3 решения, неизоморфных ~ 2
	отлично, это граф, где нет stronger petersen colouring
	это значит, что если и есть какая-то связь между peteresen colouring'ами и 2/3bm, то она не просто на уровне "все poor рёбра лежат на dominating circuits", а что-то более сложное
	всегда poor рёбра: 8-12, 10-18, 12-13, 14-19, 15-19, 16-18
	всегда rich рёбра: 0-3, 1-9, 1-5, 1-6, 2-7, 2-4, 2-9, 3-5, 3-7, 4-13, 4-5, 6-8, 6-7, 9-11
	вершины, всегда соседние с poor рёбрами: 0, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
	вершины, которые никогда не соседние с poor рёбрами: 1, 2, 3, 4, 5, 6, 7, 9
		эти вершины образуют цикл, который даже в 2bm встречается: 9 2 4 5 3 7 6 1
	интересно, что эти 2 множества не пересекаются и образуют все вершины графа

20.05g6:
	- petersen colouring: 2 решения, изоморфны друг другу вроде, то есть решение 1 неизоморфное
	это тоже граф, где нет stronger petersen colouring
	всегда poor рёбра: 8-12, 12-13
	всегда rich рёбра: 0-12, 4-13, 6-8, 8-15, 9-11, 10-13
	вершины, всегда соседние с poor рёбрами: 0
	- 2/3bm:
		если смотреть на длину dominating circuit = 19, то
		ignored бывает таким: все вершины, кроме 0, 8, 12, 13