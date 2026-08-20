# -*- coding: utf-8 -*-
"""
Puerto a dos aguas - Analisis estatico lineal con OpenSeesPy
=============================================================
CON ROTULA EN EL NODO C (momento = 0 en C)

Estructura tipo marco simetrico a dos aguas.
Apoyos articulados (pin) en A y E. Carga distribuida sobre toda la estructura.

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
    30: (6.5,  3.0),  # C_izq (lado del elem 2)
    31: (6.5,  3.0),  # C_der (lado del elem 3)
    4: (9.0,  3.0),   # D
    5: (13.0, 0.0),   # E
}

elements = {
    1: (1, 2),   # A -> B
    2: (2, 30),  # B -> C_izq
    3: (31, 4),  # C_der -> D
    4: (4, 5),   # D -> E
}

apoyos = {
    1: (1, 1, 0),  # A - articulado (pin): restringe Fx, Fy; libre Mz
    5: (1, 1, 0),  # E - articulado (pin): restringe Fx, Fy; libre Mz
}

E_acero = 2.1e6    # tf/m2
nu      = 0.30    # Poisson
G_acero = E_acero / (2.0 * (1.0 + nu))

d  = 0.320
bf = 0.300         # dimensiones areas e inercias
tf = 0.0115
tw = 0.0085

A_sec = 2.0 * bf * tf + (d - 2.0 * tf) * tw
I_sec = (bf * d**3 - (bf - tw) * (d - 2.0 * tf)**3) / 12.0

print("=" * 65)
print("  PROPIEDADES DE LA SECCION HE 340 AA")
print("=" * 65)
print(f"  A  (area)       = {A_sec*1e4:.2f} cm2")
print(f"  Ix (inercia)    = {I_sec*1e8:.1f} cm4")
print(f"  E  (modulo)     = {E_acero:.2e} tf/m2")
print()

q = 3.0  # tonf/m

# ============================================================
# 2. MODELO OpenSees
# ============================================================

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

print("=" * 65)
print("  DEFINICION DEL MODELO (CON ROTULA EN C)")
print("=" * 65)

nombres = {1: 'A', 2: 'B', 30: 'C_izq', 31: 'C_der', 4: 'D', 5: 'E'}
for tag, (x, y) in coords.items():
    ops.node(tag, x, y)
    print(f"  Nodo {tag} ({nombres[tag]}): ({x:.2f}, {y:.2f}) m")

print()

# --- ROTULA EN C ---
# Dos nodos en la misma posicion, igualar traslacion (DOFs 1,2) pero NO rotacion (DOF 3)
ops.equalDOF(30, 31, 1, 2)
print("  ROTULA en C: equalDOF(30, 31, DOFs 1,2) -> traslaciones iguales, rotaciones LIBRES")
print()

for tag, (fx, fy, mz) in apoyos.items():
    ops.fix(tag, fx, fy, mz)
    print(f"  Apoyo {tag} ({nombres[tag]}): articulado (pin)")

print()

matTag = 1
ops.uniaxialMaterial('Elastic', matTag, E_acero)

geoTag = 1
ops.geomTransf('Linear', geoTag, 0, 0, 1)

print("  Elementos:")
nombres_elem = {1: 'A-B', 2: 'B-C', 3: 'C-D', 4: 'D-E'}
for tag, (ni, nj) in elements.items():
    ops.element('elasticBeamColumn', tag, ni, nj, A_sec, E_acero, I_sec, geoTag)
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
    print(f"    Elem {tag} ({nombres_elem[tag]}): Nodos ({ni}->{nj}), L = {L:.4f} m")

print()

# ============================================================
# 3. CARGAS
# ============================================================

print("=" * 65)
print("  CARGAS APLICADAS")
print("=" * 65)

# Carga vertical global q = 3.0 tf/m hacia abajo sobre TODOS los elementos
# Para OpenSees eleLoad -beamUniform usa coord. LOCALES:
#   rotacion global->local: fx_l = cos*Fx_g + sin*Fy_g, fy_l = -sin*Fx_g + cos*Fy_g
#   Para carga (0, -q):  wx_l = -q*sin(alpha),  wy_l = -q*cos(alpha)

# Carga vertical proyectada q = 3.0 tf/m (por metro horizontal)
# En SAP2000 esto es "Projected load" en -Y global
# Para OpenSees eleLoad -beamUniform (por metro de ELEMENTO):
#   Wy_loc = -q * cos^2(alpha),  Wx_loc = -q * sin(alpha) * cos(alpha)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

for tag, (ni, nj) in elements.items():
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj - xi)**2 + (yj - yi)**2)
    cos_a = (xj - xi) / L
    sin_a = (yj - yi) / L
    Lx = xj - xi  # proyeccion horizontal
    wy_local = -q * cos_a**2           # componente local Y
    wx_local = -q * sin_a * cos_a      # componente local X
    ops.eleLoad('-ele', tag, '-type', '-beamUniform', wy_local, wx_local)
    carga_elem = q * abs(Lx)
    print(f"  Elem {tag} ({nombres_elem[tag]}): Lx={abs(Lx):.1f}m, carga={carga_elem:.1f}tf, wy_loc={wy_local:+.4f}  wx_loc={wx_local:+.4f}")

carga_total = q * 13.0  # 3.0 * (4+2.5+2.5+4) = 39.0 tf
print(f"  Carga total = q * Lx_total = {q:.1f} * 13.0 = {carga_total:.1f} tf")
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
    print("  Analisis completado exitosamente.")
print()

# ============================================================
# 5. DESPLAZAMIENTOS
# ============================================================

print("=" * 65)
print("  DESPLAZAMIENTOS EN NODOS")
print("=" * 65)
print(f"  {'Nodo':<12} {'Ux (m)':<14} {'Uy (m)':<14} {'Rot (rad)':<16}")
print("-" * 60)

tags_ord = [1, 2, 30, 31, 4, 5]
desplazamientos = {}
for tag in tags_ord:
    ux = ops.nodeDisp(tag, 1)
    uy = ops.nodeDisp(tag, 2)
    rz = ops.nodeDisp(tag, 3)
    desplazamientos[tag] = (ux, uy, rz)
    print(f"  {nombres[tag]:<12} {ux:>+.8f}   {uy:>+.8f}   {rz:>+.12f}")

# Verificar que la rotula funciona: Rot(C_izq) != Rot(C_der)
rot_cizq = desplazamientos[30][2]
rot_cder = desplazamientos[31][2]
print()
print(f"  Verificacion rotula en C:")
print(f"    Rot(C_izq) = {rot_cizq:+.10f} rad")
print(f"    Rot(C_der) = {rot_cder:+.10f} rad")
print(f"    Diferencia = {abs(rot_cizq - rot_cder):.10f} rad")
print(f"    Momento en C (debe ser ~0): Elem2_j M, Elem3_i M")

print()

# ============================================================
# 6. FUERZAS INTERNAS
# ============================================================

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

# ============================================================
# 7. REACCIONES
# ============================================================

print("=" * 65)
print("  REACCIONES EN APOYOS")
print("=" * 65)

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

print(f"  {'Punto':<8} {'Rx (tonf)':<14} {'Ry (tonf)':<14} {'Mz (tf.m)':<14}")
print("-" * 55)
for tag in apoyos:
    rx, ry, mz = reacciones[tag]
    print(f"  {nombres[tag]:<8} {rx:>+12.4f}   {ry:>+12.4f}   {mz:>+12.4f}")

print()

# ============================================================
# 8. VERIFICACION DE EQUILIBRIO
# ============================================================

Fx_total = sum(r[0] for r in reacciones.values())
Fy_total = sum(r[1] for r in reacciones.values())

print("=" * 65)
print("  VERIFICACION DE EQUILIBRIO GLOBAL")
print("=" * 65)
print(f"  SumFx = {Fx_total:+.6f} tonf  (esperado: 0)")
print(f"  SumFy = {Fy_total:+.6f} tonf  (esperado: +{carga_total:.1f})")
print()

# Equilibrio en nodos internos
print("=" * 65)
print("  VERIFICACION DE EQUILIBRIO EN NODOS INTERNOS")
print("=" * 65)
for nodo_int in [2, 4]:  # B, D
    fx_n, fy_n, mz_n = 0.0, 0.0, 0.0
    for tag, (ni, nj) in elements.items():
        forces = fuerzas_elem[tag]
        if forces is None or len(forces) < 6:
            continue
        if ni == nodo_int:
            fx_n += forces[0]; fy_n += forces[1]; mz_n += forces[2]
        elif nj == nodo_int:
            fx_n += forces[3]; fy_n += forces[4]; mz_n += forces[5]
    print(f"  {nombres[nodo_int]}: Fx={fx_n:+.6f}  Fy={fy_n:+.6f}  Mz={mz_n:+.6f}")

# En C (rotula): los momentos de cada lado deben ser ~0
m_elem2_j = fuerzas_elem[2][5] if fuerzas_elem[2] is not None else 0
m_elem3_i = fuerzas_elem[3][2] if fuerzas_elem[3] is not None else 0
print(f"  C (rotula): M_elem2_j = {m_elem2_j:+.6f}  M_elem3_i = {m_elem3_i:+.6f}  (esperado: ~0)")

print()

# ============================================================
# 9. RESUMEN
# ============================================================

print("=" * 65)
print("  RESUMEN PARA VALIDACION")
print("=" * 65)
print()
print("  Desplazamientos verticales:")
print(f"    B:  Uy = {desplazamientos[2][1]*1000:+.4f} mm")
print(f"    C:  Uy = {desplazamientos[30][1]*1000:+.4f} mm (promedio, ambos nodos)")
print(f"    D:  Uy = {desplazamientos[4][1]*1000:+.4f} mm")
print()
print(f"  Momento en C (rotula):")
print(f"    Elem 2 en j: M = {m_elem2_j:+.6f} tf.m")
print(f"    Elem 3 en i: M = {m_elem3_i:+.6f} tf.m")
print()
print("  Reacciones:")
for tag in apoyos:
    rx, ry, mz = reacciones[tag]
    print(f"    {nombres[tag]}: Ry={ry:+.4f}  Rx={rx:+.4f}  Mz={mz:+.4f}")
print(f"    Suma Ry = {reacciones[1][1]+reacciones[5][1]:+.4f} tonf")

print()
print("=" * 65)
print("  ANALISIS COMPLETADO - CON ROTULA EN C")
print("=" * 65)

# ============================================================
# 10. FUERZAS MAXIMAS POR ELEMENTO
# ============================================================

def calcular_fuerzas_al_elemento(tag, elements, fuerzas, coords, q, n_pts=100):
    """Calcula N, V, M a lo largo del elemento usando fuerzas extremo + carga."""
    ni, nj = elements[tag]
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj-xi)**2 + (yj-yi)**2)
    cos_a = (xj-xi)/L
    sin_a = (yj-yi)/L

    f = fuerzas[tag]
    Fx_i, Fy_i, Mz_i = f[0], f[1], f[2]

    Lx = xj - xi
    wy_loc = -q * cos_a**2
    w = wy_loc

    s_vals = np.linspace(0, L, n_pts)
    N_arr = np.zeros(n_pts)
    V_arr = np.zeros(n_pts)
    M_arr = np.zeros(n_pts)

    for idx, s in enumerate(s_vals):
        N_arr[idx] = Fx_i*cos_a + Fy_i*sin_a + (-q*sin_a*cos_a)*s
        V_arr[idx] = -Fx_i*sin_a + Fy_i*cos_a + w*s
        M_arr[idx] = Mz_i + (-Fx_i*sin_a + Fy_i*cos_a)*s + w*s**2/2

    return s_vals, N_arr, V_arr, M_arr

print()
print("=" * 65)
print("  FUERZAS MAXIMAS POR ELEMENTO")
print("=" * 65)

fuerzas_maximas = {}
for tag in elements:
    s_vals, N_arr, V_arr, M_arr = calcular_fuerzas_al_elemento(tag, elements, fuerzas_elem, coords, q)
    ni, nj = elements[tag]
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj-xi)**2 + (yj-yi)**2)

    result = {
        'N_max': float(np.max(N_arr)),
        'N_min': float(np.min(N_arr)),
        'V_max': float(np.max(V_arr)),
        'V_min': float(np.min(V_arr)),
        'M_max': float(np.max(M_arr)),
        'M_min': float(np.min(M_arr)),
    }
    fuerzas_maximas[nombres_elem[tag]] = result

    print(f"  Elemento {tag} ({nombres_elem[tag]}), L={L:.2f}m:")
    print(f"    Axial:   N_max={result['N_max']:+.4f} tf,  N_min={result['N_min']:+.4f} tf")
    print(f"    Cortante: V_max={result['V_max']:+.4f} tf,  V_min={result['V_min']:+.4f} tf")
    print(f"    Momento: M_max={result['M_max']:+.4f} tf.m, M_min={result['M_min']:+.4f} tf.m")
    print()

# ============================================================
# 11. EXPORTAR RESULTADOS A JSON
# ============================================================

resultados = {
    "seccion": {"tipo": "HE 340 AA", "d_mm": d*1000, "bf_mm": bf*1000,
                "tf_mm": tf*1000, "tw_mm": tw*1000, "A_cm2": round(A_sec*1e4, 2),
                "Ix_cm4": round(I_sec*1e8, 1)},
    "material": {"E_tf_m2": E_acero},
    "carga": {"q_tf_m": q, "tipo": "proyectada vertical", "carga_total_tf": q * 13.0},
    "desplazamientos_m": {nombres[k]: {"ux": v[0], "uy": v[1], "rz": v[2]}
                          for k, v in desplazamientos.items()},
    "reacciones": {nombres[k]: {"Rx": v[0], "Ry": v[1], "Mz": v[2]}
                   for k, v in reacciones.items()},
    "fuerzas_internas_extremos": {},
    "fuerzas_maximas_por_elemento": fuerzas_maximas,
}

for tag, f in fuerzas_elem.items():
    ni, nj = elements[tag]
    resultados["fuerzas_internas_extremos"][nombres_elem[tag]] = {
        "i": {"nodo": nombres[ni], "Fx": f[0], "Fy": f[1], "Mz": f[2]},
        "j": {"nodo": nombres[nj], "Fx": f[3], "Fy": f[4], "Mz": f[5]},
    }

output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados.json")
with open(output_file, "w", encoding="utf-8") as f:
    json.dump(resultados, f, indent=2, ensure_ascii=False)

print(f"  Resultados exportados a: resultados.json")
