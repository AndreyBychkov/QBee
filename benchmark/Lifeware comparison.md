# Benchmark comparison

---

| System   |      Lifeware      |  QBee | Qbee in Newton Polyhedral | QBee with Newton Polyhedral vertices |
|:----------:|:-------------:|:------:|:-------:|:-------:|
| Circular 3 |  83ms | 6ms | 13ms | 12ms |
| Circular 4 |  107ms | 105ms | 155ms | 129ms |
| Circular 5 | 570ms | 28ms | 53ms | 101ms |
| Circular 6 | 37.8s | 3.23s | 3.42s | 1.94s |
| Circular 7 | nt | 122ms | 325ms | 347ms |
| Circular 8 | nt | 41s | 42.5s | 23.4s |
| Hard 2 | ? | ? | na | 40.4s |
| Hard 3 | 1159ms* | 27.1s | na | 22.1s |
| Hard 4 | 20s* | 375.4s | na | 317s |
| Hard 5 | ? | 918s | na | 372s |
| Hill 2 | 77ms | 2ms | 3ms | 3ms |
| Hill 3 | 79ms | 2ms | 4ms | 4ms |
| Hill 4 | 81ms | 4ms | 8ms | 6ms |
| Hill 5 | 83ms | 5ms | 9ms | 8ms |
| Hill 6 | 93ms | 10ms| 19ms | 14ms |
| Hill 7 | 98ms | 13ms | 25ms | 21ms |
| Hill 8 | 127ms | 48ms | 79ms | 51ms |
| Hill 9 | 197ms | 17ms | 30ms | 37ms |
| Hill 10 | 430ms | 69ms | 104ms | 83ms |
| Hill 15 | 67s | 486ms | 831ms | 790ms |
| Hill 20 | nt | 3.36s | 8.7s | 4s |
| Long Monom 2 | 90ms* | 17ms | 29ms | 15ms |
| Long Monom 3 | 407ms* | 219s | na | 230s |

nt = not terminated: more than 2^32 iterations
na = not applicable
* order is not optimal



---

## Machine info for QBee

* **Processor**: Intel64 Family 6 Model 158 Stepping 10, GenuineIntel
* **Machine**: AMD64
* **Python**: CPython 3.8.5
* **System**: WSL Windows 10 Ubuntu 20.04


