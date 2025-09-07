# multiloops-classification

data

- raw_article_data (Dane surowe przepisane z artykułów)

laing_experiment (odtworzenie experymentu z 2012 oraz optymalizacja)

lamiable_experiment (pliki źródłowe z pracy Lamiable 2012)

graph_approach (z wykorzystaniem biblioteki forgi)

należy poprawić w bibliotece forgi następujący plik:
`.venv/lib/python3.13/site-packages/forgi/threedee/model/coarse_grain.py`
i wykomentować:

```
if key[0] in ["m", "i"]:
    # Bulge coordinates are redundant. They can be deduced from the stem coordinates.
    continue
```
