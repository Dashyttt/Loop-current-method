# Реализация алгоритма метода контурных токов на языке C++

## Цель работы:

Целью данной лабораторной работы было разработать программу для поиска контурных токов в цепи, используя матрицы.
## Задачи работы:

1. Изучить теоретические основы метода контурных токов (МКТ) и его применение в расчетах электрических цепей.
2. Разработать программу на языке программирования C++, которая реализует МКТ через матрицы.
3. Проверить корректность работы программы на нескольких тестовых примерах.
4. Сравнить результаты расчетов МКТ с ожидаемыми результатами.
В матричном виде система уравнений для МКТ выглядит следующим образом:
CZ(C)^T * I = CE, где C – матрица контуров, Z – матрица сопротивлений и E – матрица напряжений.
Алгоритм заполнения матрицы C: i–я строка соответствует независимому контуру i, а j–й столбец соответствует ветви j, причём элемент C(ij) равен
- 0, если ребро j не входит в контур i;
- 1, если ребро входит в контур, и направление ребра соответствует направлению обхода контура;
- –1, если ребро входит в контур, и направление ребра противоположно направлению обхода контура.
Матрица Z – диагональная матрица, в которой диагональный элемент Zii равен сопротивлению i–го ребра, а недиагональные элементы равны нулю. I – матрица-столбец контурных токов. E – матрица-столбец источников ЭДС, где каждый элемент равен ЭДС источника в соответствующем ребре, причём эта величина нулевая, если в данном ребре источник ЭДС отсутствует положительная, если направление ЭДС источника совпадает с направлением тока в ребре; и отрицательная в противном случае.
Решением уравнения является матрица-столбец контурных токов
I = (CZ(C)^T)^-1 * CE 

В реализации алгоритма уравнение решается методом Гаусса.