# algoritmos-t1-2026
Tarea 1 del curso Análisis de Algoritmos, Universidad de Concepción

# Compilación

```
g++ src/test.cpp -I include/ -o matrixmult
```

# Ejecución
```
./matrixmult <strassen | rowcol> <MATRIX_SIZE> [STRASSEN_THRESHOLD] [--print-result]
```

donde `STRASSEN_THRESHOLD` y `--print-result` son argumentos opcionales. El valor por defecto de `STRASSEN_THRESHOLD` es `n0 = 1`.
