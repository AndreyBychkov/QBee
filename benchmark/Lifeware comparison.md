# Benchmark comparison

---

| System   |      Lifeware      |  QBee | Qbee in Newton Polyhedral |
|:----------:|:-------------:|:------:|:-------:|
| Circular 3 |  81ms | 6ms | 11ms | 
| Circular 4 |    113ms   |   265ms | 170ms |
| Circular 5 | 878ms |    29ms | 53ms |
| Circular 6 | 54.1s |    17.68s | 3.94s |
| Hard 3 | 1560ms |    27.9s* | - |
| Hard 4 | 28.3s |    383.5s* | - |
| Hill 2 | 81ms |    2ms | 2ms |
| Hill 3| 75ms? |    2ms | 4ms |
| Hill 4 | 85ms | 4ms | 11ms |
| Hill 5 | 86ms | 5ms | 12ms |
| Hill 6 | 103ms | 11ms| 32ms |
| Hill 7 | 106ms | 13ms | 43ms |
| Hill 8 | 151ms | 47ms | 97ms |
| Hill 10 | 580ms | 58ms | 160ms
| Hill 15 | 91.6s | 0.48s | 1s |
| Long Monom 2 | 104ms | 17ms | 32ms
| Long Monom 3 | 647ms | 345s* | - 

* lower order than Lifeware(?)


---

## Machine info for QBee

* **Processor**: Intel64 Family 6 Model 158 Stepping 10, GenuineIntel
* **Machine**: AMD64
* **Python**: CPython 3.8.5
* **System**: WSL Windows 10 Ubuntu 20.04


