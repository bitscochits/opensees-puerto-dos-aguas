import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import openseespy.opensees as ops
import numpy as np
import json
import os

nombres = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}
nombres_elem = {1: 'A-B', 2: 'B-C', 3: 'C-D', 4: 'D-E'}

coords = {
    1: (0.0,  0.0),   # A
    2: (0.0,  2.0),   # B
    3: (0.0,  3.0),   # C
    4: (5.0,  2.0),   # D
    5: (8.0,  2.0),   # E
}

elements = {
    1: (1, 2),   # A -> B
    2: (2, 3),   # B -> C
    3: (3, 4),   # C -> D
    4: (4, 5),   # D -> E
}

apoyos = {
    1: (1, 1, 1),  # A - empotrado: restringe Fx, Fy, Mz
    5: (1, 1, 0),  # E - articulado (pin): restringe Fx, Fy; libre Mz
}

E_acero = 2.1e6    # tf/m2
nu      = 0.30
G_acero = E_acero / (2.0 * (1.0 + nu))

# Seccion 0.4x0.4 para columnas (A-B, B-C)
A_seca = 0.4 * 0.4
I_seca = (0.4 * 0.4**3) / 12.0

# Seccion 0.4x0.4 para vigas (C-D, D-E)
A_secb = 0.4 * 0.4
I_secb = (0.4 * 0.4**3) / 12.0

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- NODOS ---
for tag, (x, y) in coords.items():
    ops.node(tag, x, y)

# --- APOYOS ---
for tag, (fx, fy, mz) in apoyos.items():
    ops.fix(tag, fx, fy, mz)
    tipo = "empotrado" if mz == 1 else "articulado (pin)"
    print(f"  Apoyo {tag} ({nombres[tag]}): {tipo}")

print()

# --- MATERIAL Y GEOMETRIA ---
matTag = 1
ops.uniaxialMaterial('Elastic', matTag, E_acero)

geoTag = 1
ops.geomTransf('Linear', geoTag, 0, 0, 1)

# --- ELEMENTOS ---
for tag, (ni, nj) in elements.items():
    if tag <= 2:
        ops.element('elasticBeamColumn', tag, ni, nj, A_seca, E_acero, I_seca, geoTag)
    else:
        ops.element('elasticBeamColumn', tag, ni, nj, A_secb, E_acero, I_secb, geoTag)

# --- CARGAS ---
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

ops.load(4, 0.0, -20.0, 0.0)
ops.eleLoad('-ele', 1, '-type', '-beamUniform', -17.0, 0.0)
ops.eleLoad('-ele', 2, '-type', '-beamUniform', -17.0, 0.0)

ops.constraints('Plain') #manejo de restricciones (apoyos, empotramientos, etc), en este caso se usa el metodo "Plain" que es el mas simple
ops.numberer('RCM') #eficiencia en el reordenamiento de ecuaciones, elimina ceros en la matriz de rigidez
ops.system('BandGeneral') #resulve solo el sistema de ecuaciones, en este caso se usa el metodo "BandGeneral" que es el mas simple
ops.test('NormDispIncr', 1e-10, 10) #criterio de convergencia, en este caso se usa el metodo "NormDispIncr" que es el mas simple, con tolerancia 1e-10 y maximo 10 iteraciones
ops.algorithm('Linear') #lo resuelve de manera lineal, no iterativa
ops.integrator('LoadControl', 1.0) #controla la carga aplicada, en este caso se aplica toda la carga de una vez (1.0)
ops.analysis('Static')#tipo de analisis, en este caso se usa el metodo "Static" que es el mas simple
ok = ops.analyze(1)#realiza el analisis, en este caso se hace un solo paso de analisis (1)

if ok != 0:
    print("  ERROR: el analisis no convergio.")
else:
    print("  Analisis completado exitosamente.")
print()

print("=" * 65)
print("  FUERZAS INTERNAS EN ELEMENTOS (coord. globales)")
print("=" * 65)
print()

fuerzas_elem = {}
for tag in elements:
    ni, nj = elements[tag]
    forces = ops.eleForce(tag)
    fuerzas_elem[tag] = forces

    print(f"  Elemento {tag} ({nombres_elem[tag]}):")
    if forces is not None and len(forces) >= 6:
        print(f"    i (Nodo {ni}): Fx={forces[0]:>+.4f}  Fy={forces[1]:>+.4f}  Mz={forces[2]:>+.4f}")
        print(f"    j (Nodo {nj}): Fx={forces[3]:>+.4f}  Fy={forces[4]:>+.4f}  Mz={forces[5]:>+.4f}")
    print()
