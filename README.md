# Шаблон MPI-программы для параллельного решения систем линейных уравнений

## Требования

[Список задач и общие требования к решению задач](https://zenderro.github.io/programming-semester-5/presentations/LinearSystemsTasks.pdf)

## Автоматическое тестирование

Для принятия задачи необходимо, но не достаточно, чтобы программа прошла автоматические тесты. Таймаут на выполнение каждого теста составляет 1 мин.

### Оформление вывода на экран

На экран должна быть выведена исходная матрица и её правая часть в виде матрицы размера `(m + 1) x m`, где `m` - ограничение на количество выводимых символов.

Решение выводится на экран с учётом ограничения на количество выводимых символов и оформляется в отдельной строке, начинающейся с `Solution:`.

Время выполнения алгоритма должно выводиться на экран и оформляться в строке, начинающейся с `Time:`. Время алгоритма должно замеряться только для функции решения СЛУ. Расчёт погрешности и невязки не должен учитываться.

Время работы каждого процесса должно выводиться на экран и оформляться в строке, начинающейся с `Time of process i:`, где `i` - это номер процесса.

Вывод нормы невязки на экран оформляется в виде отдельной строки, начинающейся с `Residual:`.

Вывод нормы погрешности на экран оформляется в виде отдельной строки, начинающейся с `Error:`.

Тесты включают проверку возвращаемого значения. При успешном завершении программы должен возвращаться код `0`. При возникновении ошибки должен быть возвращён один из кодов ошибок.

### Коды ошибок

`-1`. Некорректные аргументы:

   - неправильное количество аргументов;
   - некорректные значения аргументов (`n<=0`, `m<=0`, `m > n`)
   - не удалось считать число для `n`, `m` или `k`;
   - некорректное значение `k`;
   - при `k == 0` не указан файл.

`-2`. Не удалось выделить память.

`-3`. Ошибки чтения файла:

   - указанного файла не существует;
   - файл пуст;
   - не удалось считать матрицу из файла.

`-4`. Вырожденная матрица.

## Общая информация по шаблону

Заготовка программы находится в файле `main.c`.
Проще всего писать программу в этом файле. При необходимости можно его переименовать.

Программа может состоять из нескольких модулей.
При проверке компилируются все файлы с расширением `.c` в этом каталоге.

Поэтому перед отправкой нужно убедиться, что в последнем наборе изменений
записаны только нужные файлы

Программу можно писать на языке C++, если требования к заданию это допускают. Для этого нужно переименовать файл `main.c` в `main.cpp`. Для программ на языке C++ в репозитории должны быть только файлы исходного кода с расширением `.cpp`. Для программ на языке C — только файлы исходного кода с расширением `.c`.

## Особенности MPI

Шпаргалка с основными командами для передачи данных между процессами с использованием MPI доступна по [этой ссылке](https://zenderro.github.io/programming-semester-6/MPI-cheatsheet).

### Дополнительные требования

В силу того, что технология MPI предназначена для работы с системами с распределённой памятью, выдвигаются дополнительные требования к программам.

1. Не более O(n) пересылок, объём каждой — не более O(n). Это означает, что всю матрицу целиком пересылать внутри метода нельзя.
2. Матрица должна быть распределена между процессами. Если процессов больше одного, ни один процесс не должен выделять память под всю матрицу (или под всю обратную). Использование памяти для P процессов должно быть n^2/P+O(n) для решения линейной системы и 2n^2/P+O(n) для обращения матрицы.
3. Чтение матрицы из файла и генерация матрицы должно осуществляться одним процессом.
4. Вывод матрицы на экран должно осуществляться одним процессом.

## Сборка

Сборка программы осуществляется с использованием `Makefile`.

Предусмотрено три режима сборки:

- `make` - "быстрая" сборка с оптимизацией, используется для измерения времени работы программы на больших размерах матрицы;
- `make debug` - отладочная сборка для отладки с использованием `gdb`;
- `make test` - сборка для тестирования с проверками на на утечку памяти и на выход за границы массива.  Эта версия программы будет использоваться запуска тестов.

Для перехода между разными режимами необходимо выполнить команду `make clean`, которая удаляет все файлы, созданные во время сборки.

## Проверка

В дополнение к автотестам в составе задания, шаблон содержит следующие общие проверки:

- `build` — собирает программу
    с проверками на утечку памяти и на выход за границы массива.
    Эта версия программы будет использоваться запуска тестов.
    Можно запустить на своём компьютере (локально) командой `make test`

    В Ubuntu установить нужные программы, включая `git`, можно командой:

```
sudo apt-get install -y build-essential
```

- `lint` — выполняет стандартные проверки:
    - компиляция со всеми предупреждениями
    - проверка кода анализаторами `clang-tidy` и `cppcheck`
    - проверка форматирования кода анализатором `clang-format` (нужна версия 10)
    Если перечисленные средства установлены,
    проверку можно запустить локально командой `.github/lint.sh`.

    В Ubuntu установить нужные программы можно командой:

```
sudo apt-get install -y clang-tidy clang-format cppcheck
```

Некоторые замечания по результатам проверки `clang-tidy` можно исправить автоматически.
Для этого вызовите команду:

```
.github/clang-tidy.sh --fix-errors main.c
```

Файлы в каталоге `.github` студентам редактировать запрещено.

## Требования к форматированию кода

Не все требования к оформлению можно проверить автоматически.
Даже к программе, которая проходит проверку `clang-format` преподаватель может оставить замечания по оформлению, стилю идентификаторов и удобству чтения.
Общий стиль оформления близок к K&R.

Размер отступов — 4 пробела. Внутри каждого блока (внутри фигурных скобок) добавляется новый уровень отступа.

Фигурные скобки, открывающие блоки, пишутся на той же строке, что и операторы `if`, `for`, `while`, `do`, `switch`, через пробел.
На отдельной строке фигурные скобки пишутся в определении функций и пространств имён в C++.
Перед скобками после операторов `if`, `for`, `while`, `switch` ставится пробел.
Перед скобками с аргументами в вызове функции пробел НЕ ставится.

Если `clang-format` есть только более старых версий,
можно автоматически переформатировать код следующей командой.
Внимание: ключ `-i` переписывает файлы. Это может «испортить» форматирование, если в программе есть синтаксические ошибки.
Поэтому на всякий случай, сделайте коммит перед запуском `clang-format`.

```
.github/clang-format.sh -i ./main.c
```

Пример правильно отформатированного кода:

```c
#include <stdio.h>

void output(char character, char *line);

void output(char character, char *line)
{
    printf("Argument: %s\nFirst char: %c", line, character);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Not enough arguments: %d!\n", argc);
        return 1;
    }

    output(argv[1], argv[1][0]);
    return 0;
}
```
