# Benchmark comparison

---

| System   |      Lifeware      |  QBee | Qbee in Newton Polyhedral | QBee with Newton Polyhedral vertices |
|:----------:|:-------------:|:------:|:-------:|:-------:|
| Circular 3 |  83ms | 6ms | 13ms | 12ms |
| Circular 4 |  107ms | 261ms | 170ms | 383ms |
| Circular 5 | 570ms | 28ms | 53ms | 435ms |
| Circular 6 | 37.8s | 17.1s | 3.94s | 12.7s |
| Circular 7 | nt | 220ms | 325ms | 3.94s |
| Hard 3 | 1159ms* | 27.1s | na | 24.1s |
| Hard 4 | 20s* | 375.4s | na | 356s |
| Hill 2 | 77ms | 2ms | 3ms | 3ms |
| Hill 3 | 79ms | 3ms | 4ms | 4ms |
| Hill 4 | 81ms | 4ms | 8ms | 6ms |
| Hill 5 | 83ms | 5ms | 9ms | 8ms |
| Hill 6 | 93ms | 11ms| 19ms | 4ms |
| Hill 7 | 98ms | 14ms | 25ms | 21ms |
| Hill 8 | 127ms | 48ms | 79ms | 59ms |
| Hill 9 | 197ms | 18ms | 30ms | 47ms |
| Hill 10 | 430ms | 69ms | 104ms | 84ms |
| Hill 15 | 67s | 486ms | 831ms | 826ms |
| Hill 20 | nt | 3.36s | 8.7s | 4.2s |
| Long Monom 2 | 90ms* | 18ms | 29ms | 15ms |
| Long Monom 3 | 407ms* | 345s | na | 335.6s |

nt = not terminated: more than 2^32 iterations
na = not applicable
* order is not optimal



---

## Machine info for QBee

* **Processor**: Intel64 Family 6 Model 158 Stepping 10, GenuineIntel
* **Machine**: AMD64
* **Python**: CPython 3.8.5
* **System**: WSL Windows 10 Ubuntu 20.04


