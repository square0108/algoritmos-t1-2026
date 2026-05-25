# algoritmos-t1-2026
El presente repositorio contiene una implementación básica del algoritmo de multiplicación de matrices de Strassen, ofreciendo un threshold ajustable de tamaño de matriz donde el programa deja de ejecutar Strassen y procede a ejecutar multiplicación fila por columna.

# Compilación

```
g++ -std=c++17 src/test.cpp -I include/ -o matrixmult
```

# Ejecución
```
./matrixmult <strassen | rowcol> <MATRIX_SIZE> [STRASSEN_THRESHOLD] [--print-result]
```

donde `STRASSEN_THRESHOLD` y `--print-result` son argumentos opcionales. El valor por defecto de `STRASSEN_THRESHOLD` es `n0 = 1`.
