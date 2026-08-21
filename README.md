# Analisis Estructural con OpenSeesPy

Proyecto de analisis estatico lineal de estructuras usando OpenSeesPy (Python 3.10).

## Requisitos

- Python 3.10 (no compatible con 3.14)
- OpenSeesPy 3.5.1.3
- NumPy

## Instalacion

```bash
# Crear entorno virtual
py -3.10 -m venv .venv

# Activar entorno
# Windows:
.venv\Scripts\activate
# Linux/Mac:
source .venv/bin/activate

# Instalar dependencias
pip install openseespy numpy
```

## Archivos

| Archivo | Descripcion |
|---------|-------------|
| `frame_analysis_rotula.py` | Marco a dos aguas con rotula en C, apoyos articulados, carga proyectada |
| `ejercicio_ayudantia.py` | Ejercicio de ayudantia - marco con empotramiento parcial |

## Ejecucion

```bash
py -3.10 frame_analysis_rotula.py
py -3.10 ejercicio_ayudantia.py
```

## Unidades

- Fuerza: tonf (tonelada-fuerza)
- Longitud: m (metro)
