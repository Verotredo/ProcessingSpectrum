## Выбор метода классификации
На вход метода подаются большое количество наблюдений со 120 признаками, у которых заранее могут быть известны классы, к которому они принадлежат.
В связи с биологическим характером задачи, признаки не являются независимыми. 
В результате проведения классификации наблюдения будут приписаны к одному из классов - родов цианнобактерий. 
<br>
Для данной задачи требованиями к методам являются:
+  устойчивость результатов классификации 
(при повторном проведении классификации при тех же параметрах были бы получены те же результаты); 
+  применение в методе машинного обучения;
+  простота реализации;
+  работа с априорной информацией;
+  использование статистических оценок;
+  Минимальное влияние от требования к независимости признаков, к распределениям признаков;
+  Работа метода с категориальными признаками;
+  Возможность визуального представления результата
+  Работа с большим объемом выборок без особых модификации методов.
<br>
Рассмотрим различия методов в таблице. <br>
| Метод                            | Устойчивость | МО   | Простота | Априорная инф-ция | Стат. оценка | 
| -------------------------------- | :----------: | :----: | :------: | :---------------: | :----------: | 
| Метод k-средних                  | - | - |    +     |       -           |      -       |
| Иерархическая кластеризация      | + | - |    +     |       -           |      -       |  
| Деревья классификации            | + | + |    -     |       -           |      -       |            
| Метод k ближайших соседей        | - | - |    +     |       +           |      -       |
| Наивный байесовский классификатор| - | + |    -     |       +           |      +       |
| Линейный дискриминантный анализ  |      +       | - |    -     |       +           |      +       |

<br>
| Метод                            | Без требований к пр-кам | Катег. пр-ки | Визуал. представление | Большие объемы |
| -------------------------------- | :---------------------: | :----------: | :-------------------: | :------------: |
| Метод k-средних                  |			+			 |    -         |          +            |       -        | 
| Иерархическая кластеризация      |			+			 |    +         |          +            |       +        |  
| Деревья классификации            |			+			 |    +         |          +            |       +        |             
| Метод k ближайших соседей        |			+			 |    -         |          +            |       -        |
| Наивный байесовский классификатор|			-			 |    -         |          -            |       +        |
| Линейный дискриминантный анализ  |			*			 |    +         |          *            |       +        |

<br>
По результатам сравнения был выбран линейный дискриминантный анализ, т.к. его результаты классификации стабильны, в отличии от байесовского классификатора, посредством канонического анализа результаты
можно представить на графике, а также он использует априорную информацию - наблюдаемые классы и использует в оценке статистические методы.

## Источники (потом перепишем)
1. jain2010.pdf

2. Ref4_k-means.pdf

3. johnson1967.pdf

