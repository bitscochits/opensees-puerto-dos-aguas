# -*- coding: utf-8 -*-
"""
Puerto a dos aguas - Analisis estatico lineal con OpenSeesPy
=============================================================
Estructura tipo marco simetrico a dos aguas.
Apoyos empotrados en A y E. Carga distribuida sobre la viga superior B-C-D.

Unidades: tonf (fuerza), m (longitud)

Seccion HE 340 AA:
  d=320mm, bf=300mm, tf=11.5mm, tw=8.5mm
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import openseespy.opensees as ops
import numpy as np
import json
import os

# ============================================================
# 1. DATOS DE ENTRADA
# ============================================================

coords = {
    1: (0.0,  0.0),   # A
    2: (4.0,  3.0),   # B
    3: (6.5,  3.0),   # C
    4: (9.0,  3.0),   # D
    5: (13.0, 0.0),   # E
}

elements = {
    1: (1, 2),  # A -> B
    2: (2, 3),  # B -> C
    3: (3, 4),  # C -> D
    4: (4, 5),  # D -> E
}

apoyos = {
    1: (1, 1, 1),  # A - empotrado
    5: (1, 1, 1),  # E - empotrado
}

E_acero = 2.1e6    # tf/m2
nu      = 0.30
G_acero = E_acero / (2.0 * (1.0 + nu))

d  = 0.320
bf = 0.300
tf = 0.0115
tw = 0.0085

A_sec = 2.0 * bf * tf + (d - 2.0 * tf) * tw
I_sec = (bf * d**3 - (bf - tw) * (d - 2.0 * tf)**3) / 12.0
Iyy_sec = (2.0 * tf * bf**3 + (d - 2.0 * tf) * tw**3) / 12.0
J_sec = (2.0 * bf * tf**3 + (d - 2.0 * tf) * tw**3) / 3.0

print("=" * 65)
print("  PROPIEDADES DE LA SECCION HE 340 AA")
print("=" * 65)
print(f"  d  (altura)     = {d*1000:.0f} mm")
print(f"  bf (ancho patin)= {bf*1000:.0f} mm")
print(f"  tf (esp. patin) = {tf*1000:.1f} mm")
print(f"  tw (esp. alma)  = {tw*1000:.1f} mm")
print(f"  A  (area)       = {A_sec*1e4:.2f} cm2")
print(f"  Ix (inercia)    = {I_sec*1e8:.1f} cm4")
print(f"  Iy (inercia yy) = {Iyy_sec*1e8:.1f} cm4")
print(f"  J  (torsion)    = {J_sec*1e8:.2f} cm4")
print(f"  E  (modulo)     = {E_acero:.2e} tf/m2")
print(f"  G  (cortante)   = {G_acero:.2e} tf/m2")
print()

q = 3.0  # tonf/m

# ============================================================
# 2. MODELO OpenSees
# ============================================================

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

print("=" * 65)
print("  DEFINICION DEL MODELO")
print("=" * 65)

nombres = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E'}
for tag, (x, y) in coords.items():
    ops.node(tag, x, y)
    print(f"  Nodo {tag} ({nombres[tag]}): ({x:.2f}, {y:.2f}) m")

print()

for tag, (fx, fy, mz) in apoyos.items():
    ops.fix(tag, fx, fy, mz)
    print(f"  Apoyo {tag} ({nombres[tag]}): empotrado [Fx=RE, Fy=RE, Mz=RE]")

print()

matTag = 1
ops.uniaxialMaterial('Elastic', matTag, E_acero)
print(f"  Material #{matTag}: Acero E = {E_acero:.2e} tf/m2")

geoTag = 1
ops.geomTransf('Linear', geoTag, 0, 0, 1)
print(f"  Transformacion geometrica #{geoTag}: Linear (plano XY)")

print()
print("  Elementos:")
nombres_elem = {1: 'A-B', 2: 'B-C', 3: 'C-D', 4: 'D-E'}
for tag, (ni, nj) in elements.items():
    ops.element('elasticBeamColumn', tag, ni, nj, A_sec, E_acero, I_sec, geoTag)
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
    print(f"    Elem {tag} ({nombres_elem[tag]}): L = {L:.4f} m")

print()

# ============================================================
# 3. CARGAS
# ============================================================

print("=" * 65)
print("  CARGAS APLICADAS")
print("=" * 65)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.eleLoad('-ele', 2, '-type', '-beamUniform', -q, 0.0)
ops.eleLoad('-ele', 3, '-type', '-beamUniform', -q, 0.0)

print(f"  Carga distribuida q = {q:.1f} tonf/m (hacia abajo, -Y)")
print(f"    Aplicada sobre elemento 2 (B-C) y elemento 3 (C-D)")
print(f"    Longitud cargada total = 5.00 m")
print(f"    Carga total = {q * 5.0:.1f} tonf")
print()

# ============================================================
# 4. ANALISIS
# ============================================================

print("=" * 65)
print("  ANALISIS")
print("=" * 65)

ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1e-10, 10)
ops.algorithm('Linear')
ops.integrator('LoadControl', 1.0)
ops.analysis('Static')
ok = ops.analyze(1)

if ok != 0:
    print("  ERROR: el analisis no convergio.")
else:
    print("  Analisis estatico lineal completado exitosamente.")
print()

# ============================================================
# 5. DESPLAZAMIENTOS
# ============================================================

print("=" * 65)
print("  DESPLAZAMIENTOS EN NODOS")
print("=" * 65)
print(f"  {'Nodo':<8} {'Punto':<8} {'Ux (m)':<14} {'Uy (m)':<14} {'Rot (rad)':<14}")
print("-" * 65)

desplazamientos = {}
for tag in coords:
    ux = ops.nodeDisp(tag, 1)
    uy = ops.nodeDisp(tag, 2)
    rz = ops.nodeDisp(tag, 3)
    desplazamientos[tag] = (ux, uy, rz)
    print(f"  {tag:<8} {nombres[tag]:<8} {ux:>+.8f}   {uy:>+.8f}   {rz:>+.10f}")

print()

# ============================================================
# 6. FUERZAS INTERNAS EN ELEMENTOS
# ============================================================

print("=" * 65)
print("  FUERZAS INTERNAS EN ELEMENTOS (coord. globales)")
print("=" * 65)
print()

# Almacenar fuerzas por nodo para calcular reacciones
fuerzas_elem = {}

for tag in elements:
    ni, nj = elements[tag]
    forces = ops.eleForce(tag)
    fuerzas_elem[tag] = forces

    print(f"  Elemento {tag} ({nombres_elem[tag]}):")
    if forces is not None and len(forces) >= 6:
        print(f"    i (Nodo {ni}): N={forces[0]:>+.4f} tf  V={forces[1]:>+.4f} tf  M={forces[2]:>+.4f} tf.m")
        print(f"    j (Nodo {nj}): N={forces[3]:>+.4f} tf  V={forces[4]:>+.4f} tf  M={forces[5]:>+.4f} tf.m")
    else:
        print(f"    No se pudieron obtener las fuerzas.")
    print()

# ============================================================
# 7. REACCIONES (calculadas desde fuerzas internas)
# ============================================================

print("=" * 65)
print("  REACCIONES EN APOYOS (calculadas desde eleForce)")
print("=" * 65)

# eleForce para elasticBeamColumn en OpenSeesPy retorna fuerzas en
# coordenadas GLOBALES: [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]
# La reaccion en un apoyo = suma de eleForce en los extremos conectados

reacciones = {}
for nodo_apoyo in apoyos:
    fx_sum = 0.0
    fy_sum = 0.0
    mz_sum = 0.0
    for tag, (ni, nj) in elements.items():
        forces = fuerzas_elem[tag]
        if forces is None or len(forces) < 6:
            continue
        if ni == nodo_apoyo:
            fx_sum += forces[0]
            fy_sum += forces[1]
            mz_sum += forces[2]
        elif nj == nodo_apoyo:
            fx_sum += forces[3]
            fy_sum += forces[4]
            mz_sum += forces[5]
    reacciones[nodo_apoyo] = (fx_sum, fy_sum, mz_sum)

print(f"  {'Nodo':<8} {'Punto':<8} {'Rx (tonf)':<14} {'Ry (tonf)':<14} {'Mz (tf.m)':<14}")
print("-" * 65)
for tag in apoyos:
    rx, ry, mz = reacciones[tag]
    print(f"  {tag:<8} {nombres[tag]:<8} {rx:>+12.4f}   {ry:>+12.4f}   {mz:>+12.4f}")

print()

# ============================================================
# 8. VERIFICACION DE EQUILIBRIO
# ============================================================

Fx_total = sum(r[0] for r in reacciones.values())
Fy_total = sum(r[1] for r in reacciones.values())
M_total  = sum(r[2] for r in reacciones.values())
carga_total = q * 5.0

print("=" * 65)
print("  VERIFICACION DE EQUILIBRIO GLOBAL")
print("=" * 65)
print(f"  SumFx (reacciones) = {Fx_total:+.6f} tonf  (esperado: 0)")
print(f"  SumFy (reacciones) = {Fy_total:+.6f} tonf  (esperado: +{carga_total:.1f} tonf)")
print(f"  Carga total apl.   = {carga_total:.1f} tonf")
print(f"  Diferencia Fy      = {abs(Fy_total - carga_total):.6f} tonf")
print()

# Verificar equilibrio en nodos internos
print("=" * 65)
print("  VERIFICACION DE EQUILIBRIO EN NODOS INTERNOS")
print("=" * 65)
for nodo_int in [2, 3, 4]:  # B, C, D
    fx_n = 0.0
    fy_n = 0.0
    mz_n = 0.0
    for tag, (ni, nj) in elements.items():
        forces = fuerzas_elem[tag]
        if forces is None or len(forces) < 6:
            continue
        if ni == nodo_int:
            fx_n += forces[0]
            fy_n += forces[1]
            mz_n += forces[2]
        elif nj == nodo_int:
            fx_n += forces[3]
            fy_n += forces[4]
            mz_n += forces[5]
    print(f"  Nodo {nodo_int} ({nombres[nodo_int]}): SumFx={fx_n:+.6f}  SumFy={fy_n:+.6f}  SumMz={mz_n:+.6f}  (esperado: 0, 0, 0)")

print()

# ============================================================
# 9. RESUMEN PARA VALIDACION
# ============================================================

print("=" * 65)
print("  RESUMEN PARA VALIDACION")
print("=" * 65)
print()
print("  Desplazamientos verticales:")
print(f"    Nodo B (cumbrera izq): Uy = {desplazamientos[2][1]*1000:+.4f} mm")
print(f"    Nodo C (cumbrera ctr): Uy = {desplazamientos[3][1]*1000:+.4f} mm")
print(f"    Nodo D (cumbrera der): Uy = {desplazamientos[4][1]*1000:+.4f} mm")
print()
print("  Reacciones verticales:")
for tag in apoyos:
    rx, ry, mz = reacciones[tag]
    print(f"    {nombres[tag]}: Ry = {ry:+.4f} tonf,  Rx = {rx:+.4f} tonf,  Mz = {mz:+.4f} tf.m")
ry_a = reacciones[1][1]
ry_e = reacciones[5][1]
print(f"    Suma Ry = {ry_a + ry_e:+.4f} tonf (carga total = {carga_total:.1f} tonf)")
print()

# ============================================================
# 10. EXPORTAR RESULTADOS A JSON
# ============================================================

resultados = {
    "seccion": {
        "tipo": "HE 340 AA",
        "d_mm": d*1000, "bf_mm": bf*1000,
        "tf_mm": tf*1000, "tw_mm": tw*1000,
        "A_cm2": round(A_sec*1e4, 2),
        "Ix_cm4": round(I_sec*1e8, 1),
    },
    "material": {
        "E_tf_m2": E_acero,
        "G_tf_m2": round(G_acero, 0),
    },
    "carga": {
        "q_tf_m": q,
        "longitud_cargada_m": 5.0,
        "carga_total_tf": q * 5.0,
    },
    "desplazamientos_m": {nombres[k]: {"ux": v[0], "uy": v[1], "rz": v[2]} for k, v in desplazamientos.items()},
    "reacciones": {nombres[k]: {"Rx": v[0], "Ry": v[1], "Mz": v[2]} for k, v in reacciones.items()},
}

output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados.json")
with open(output_file, "w", encoding="utf-8") as f:
    json.dump(resultados, f, indent=2, ensure_ascii=False)

print(f"  Resultados exportados a: resultados.json")
print()
print("=" * 65)
print("  ANALISIS COMPLETADO")
print("=" * 65)
