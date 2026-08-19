# -*- coding: utf-8 -*-
"""
Diagramas de esfuerzos y deformaciones - Puerto a dos aguas
Genera PDF con: Momento flector, Cortante, Fuerza axial, Deformada
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import openseespy.opensees as ops
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch

# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

def run_analysis(q_load=3.0, con_rotula=True):
    """Ejecuta el analisis y retorna resultados."""
    coords = {
        1: (0.0, 0.0), 2: (4.0, 3.0), 30: (6.5, 3.0),
        31: (6.5, 3.0), 4: (9.0, 3.0), 5: (13.0, 0.0),
    }
    elements = {1: (1, 2), 2: (2, 30), 3: (31, 4), 4: (4, 5)}
    apoyos = {1: (1, 1, 0), 5: (1, 1, 0)}
    q = q_load

    E_acero = 2.1e6
    d, bf, tf_s, tw = 0.320, 0.300, 0.0115, 0.0085
    A_sec = 2*bf*tf_s + (d-2*tf_s)*tw
    I_sec = (bf*d**3 - (bf-tw)*(d-2*tf_s)**3)/12.0

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    for tag, (x, y) in coords.items():
        ops.node(tag, x, y)

    if con_rotula:
        ops.equalDOF(30, 31, 1, 2)

    for tag, (fx, fy, mz) in apoyos.items():
        ops.fix(tag, fx, fy, mz)

    ops.uniaxialMaterial('Elastic', 1, E_acero)
    ops.geomTransf('Linear', 1, 0, 0, 1)

    for tag, (ni, nj) in elements.items():
        ops.element('elasticBeamColumn', tag, ni, nj, A_sec, E_acero, I_sec, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)

    for tag, (ni, nj) in elements.items():
        xi, yi = coords[ni]
        xj, yj = coords[nj]
        L = np.sqrt((xj-xi)**2 + (yj-yi)**2)
        cos_a = (xj-xi)/L
        sin_a = (yj-yi)/L
        Lx = xj - xi
        wy_loc = -q * cos_a**2
        wx_loc = -q * sin_a * cos_a
        ops.eleLoad('-ele', tag, '-type', '-beamUniform', wy_loc, wx_loc)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-10, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    ops.analyze(1)

    nodos = {}
    for tag in coords:
        nodos[tag] = {
            'ux': ops.nodeDisp(tag, 1),
            'uy': ops.nodeDisp(tag, 2),
            'rz': ops.nodeDisp(tag, 3),
        }

    fuerzas = {}
    for tag in elements:
        fuerzas[tag] = ops.eleForce(tag)

    return coords, elements, nodos, fuerzas, q

def get_element_forces_at_x(tag, elements, fuerzas, coords, q, n_pts=50):
    """Obtiene N, V, M a lo largo del elemento usando las fuerzas extremo y la carga."""
    ni, nj = elements[tag]
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    L = np.sqrt((xj-xi)**2 + (yj-yi)**2)
    cos_a = (xj-xi)/L
    sin_a = (yj-yi)/L

    f = fuerzas[tag]
    Fx_i, Fy_i, Mz_i = f[0], f[1], f[2]
    Fx_j, Fy_j, Mz_j = f[3], f[4], f[5]

    # Carga local
    Lx = xj - xi
    wy_loc = -q * cos_a**2
    wx_loc = -q * sin_a * cos_a
    w = wy_loc  # carga transversal local

    s_vals = np.linspace(0, L, n_pts)
    N_arr = np.zeros(n_pts)
    V_arr = np.zeros(n_pts)
    M_arr = np.zeros(n_pts)

    for idx, s in enumerate(s_vals):
        N_arr[idx] = Fx_i*cos_a + Fy_i*sin_a + wx_loc*s
        V_arr[idx] = -Fx_i*sin_a + Fy_i*cos_a + w*s
        M_arr[idx] = Mz_i + (-Fx_i*sin_a + Fy_i*cos_a)*s + w*s**2/2

    # Coordenadas globales a lo largo del elemento
    x_glob = xi + (xj-xi)*s_vals/L
    y_glob = yi + (yj-yi)*s_vals/L

    return s_vals, N_arr, V_arr, M_arr, x_glob, y_glob, L

# ============================================================
# EJECUTAR ANALISIS
# ============================================================

print("Ejecutando analisis...")
coords, elements, nodos, fuerzas, q = run_analysis(q_load=3.0, con_rotula=True)
print("Analisis completado.")

nombres_elem = {1: 'A-B', 2: 'B-C', 3: 'C-D', 4: 'D-E'}
nombres_nodo = {1: 'A', 2: 'B', 30: 'C1', 31: 'C2', 4: 'D', 5: 'E'}

# ============================================================
# GENERAR PDF
# ============================================================

print("Generando PDF...")

with PdfPages('diagramas_puerto.pdf') as pdf:

    # ========================================
    # 1. ESTRUCTURA ORIGINAL + DEFORMADA
    # ========================================
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_aspect('equal')
    ax.set_title('Estructura Original y Deformada\n(Puerto a dos aguas - Rotula en C, q=3 tf/m proyectada)', fontsize=13, fontweight='bold')

    esc_def = 20  # factor de escala para deformada

    # Estructura original (lineas grises punteadas)
    tags_ord = [1, 2, 30, 31, 4, 5]
    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        if t1 == 30 and t2 == 31:
            continue  # saltar la union de rotula
        x1, y1 = coords[t1]
        x2, y2 = coords[t2]
        ax.plot([x1, x2], [y1, y2], 'k--', linewidth=1, alpha=0.4, label='Original' if tag_i==0 else '')

    # Nodos originales
    for tag in [1, 2, 30, 4, 5]:
        x, y = coords[tag]
        ax.plot(x, y, 'ko', markersize=5)

    # Rotula en C
    ax.plot(coords[30][0], coords[30][1], 'o', color='red', markersize=10, markerfacecolor='white', markeredgewidth=2, label='Rotula C')

    # Estructura deformada
    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        x1d = coords[t1][0] + nodos[t1]['ux']*esc_def
        y1d = coords[t1][1] + nodos[t1]['uy']*esc_def
        x2d = coords[t2][0] + nodos[t2]['ux']*esc_def
        y2d = coords[t2][1] + nodos[t2]['uy']*esc_def
        ax.plot([x1d, x2d], [y1d, y2d], 'b-', linewidth=2, label='Deformada (x20)' if tag_i==0 else '')

    # Apoyos
    for tag in [1, 5]:
        x, y = coords[tag]
        ax.plot(x, y, 's', color='green', markersize=12, markeredgewidth=2, label='Apoyo articulado' if tag==1 else '')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.text(6.5, -1.5, f'Deformada amplificada x{esc_def}', ha='center', fontsize=10, style='italic')

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

    # ========================================
    # 2. DIAGRAMA DE MOMENTO FLECTOR
    # ========================================
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_aspect('equal')
    ax.set_title('Diagrama de Momento Flector Mz\n(Puerto a dos aguas - Rotula en C, q=3 tf/m)', fontsize=13, fontweight='bold')

    esc_m = 0.15  # escala para el diagrama de momento

    # Dibujar estructura original
    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        if t1 == 30 and t2 == 31:
            continue
        x1, y1 = coords[t1]
        x2, y2 = coords[t2]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

    max_moment = 0

    for tag in elements:
        s_vals, N_arr, V_arr, M_arr, x_glob, y_glob, L = get_element_forces_at_x(tag, elements, fuerzas, coords, q)
        ni, nj = elements[tag]
        xi, yi = coords[ni]
        xj, yj = coords[nj]
        dx, dy = xj-xi, yj-yi
        L_elem = np.sqrt(dx**2+dy**2)
        nx, ny = -dy/L_elem, dx/L_elem  # normal perpendicular al elemento

        max_moment = max(max_moment, np.max(np.abs(M_arr)))

        # Linea base
        ax.plot([xi, xj], [yi, yj], 'k-', linewidth=0.5)

        # Diagrama
        x_base = np.linspace(xi, xj, len(M_arr))
        y_base = np.linspace(yi, yj, len(M_arr))
        x_diag = x_base + M_arr * esc_m * nx
        y_diag = y_base + M_arr * esc_m * ny

        ax.fill(x_diag, y_diag, alpha=0.3, color='blue')
        ax.plot(x_diag, y_diag, 'b-', linewidth=1.5)

        # Etiquetas en extremos
        ax.text(xi + M_arr[0]*esc_m*nx*1.3, yi + M_arr[0]*esc_m*ny*1.3,
                f'{M_arr[0]:+.2f}', fontsize=8, ha='center', va='center', color='blue', fontweight='bold')
        ax.text(xj + M_arr[-1]*esc_m*nx*1.3, yj + M_arr[-1]*esc_m*ny*1.3,
                f'{M_arr[-1]:+.2f}', fontsize=8, ha='center', va='center', color='blue', fontweight='bold')

    # Etiquetas de nodos
    for tag in [1, 2, 30, 4, 5]:
        ax.annotate(nombres_nodo[tag], (coords[tag][0], coords[tag][1]),
                    textcoords="offset points", xytext=(0, 12), ha='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.text(0.02, 0.02, 'Momento en tf.m (+ = traccion lado izq/direcc. normal)\nValores en extremos de cada elemento',
            transform=ax.transAxes, fontsize=9, va='bottom', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

    # ========================================
    # 3. DIAGRAMA DE CORTANTE
    # ========================================
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_aspect('equal')
    ax.set_title('Diagrama de Fuerza Cortante V\n(Puerto a dos aguas - Rotula en C, q=3 tf/m)', fontsize=13, fontweight='bold')

    esc_v = 0.3

    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        if t1 == 30 and t2 == 31:
            continue
        x1, y1 = coords[t1]
        x2, y2 = coords[t2]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

    for tag in elements:
        s_vals, N_arr, V_arr, M_arr, x_glob, y_glob, L = get_element_forces_at_x(tag, elements, fuerzas, coords, q)
        ni, nj = elements[tag]
        xi, yi = coords[ni]
        xj, yj = coords[nj]
        dx, dy = xj-xi, yj-yi
        L_elem = np.sqrt(dx**2+dy**2)
        nx, ny = -dy/L_elem, dx/L_elem

        x_base = np.linspace(xi, xj, len(V_arr))
        y_base = np.linspace(yi, yj, len(V_arr))
        x_diag = x_base + V_arr * esc_v * nx
        y_diag = y_base + V_arr * esc_v * ny

        ax.fill(x_diag, y_diag, alpha=0.3, color='red')
        ax.plot(x_diag, y_diag, 'r-', linewidth=1.5)

        ax.text(xi + V_arr[0]*esc_v*nx*1.3, yi + V_arr[0]*esc_v*ny*1.3,
                f'{V_arr[0]:+.2f}', fontsize=8, ha='center', va='center', color='red', fontweight='bold')
        ax.text(xj + V_arr[-1]*esc_v*nx*1.3, yj + V_arr[-1]*esc_v*ny*1.3,
                f'{V_arr[-1]:+.2f}', fontsize=8, ha='center', va='center', color='red', fontweight='bold')

    for tag in [1, 2, 30, 4, 5]:
        ax.annotate(nombres_nodo[tag], (coords[tag][0], coords[tag][1]),
                    textcoords="offset points", xytext=(0, 12), ha='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.text(0.02, 0.02, 'Cortante en tf (en coord. locales)\nValores en extremos de cada elemento',
            transform=ax.transAxes, fontsize=9, va='bottom', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

    # ========================================
    # 4. DIAGRAMA DE FUERZA AXIAL
    # ========================================
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_aspect('equal')
    ax.set_title('Diagrama de Fuerza Axial N\n(Puerto a dos aguas - Rotula en C, q=3 tf/m)', fontsize=13, fontweight='bold')

    esc_n = 0.2

    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        if t1 == 30 and t2 == 31:
            continue
        x1, y1 = coords[t1]
        x2, y2 = coords[t2]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1.5)

    for tag in elements:
        s_vals, N_arr, V_arr, M_arr, x_glob, y_glob, L = get_element_forces_at_x(tag, elements, fuerzas, coords, q)
        ni, nj = elements[tag]
        xi, yi = coords[ni]
        xj, yj = coords[nj]
        dx, dy = xj-xi, yj-yi
        L_elem = np.sqrt(dx**2+dy**2)
        tx, ty = dx/L_elem, dy/L_elem  # tangente

        x_base = np.linspace(xi, xj, len(N_arr))
        y_base = np.linspace(yi, yj, len(N_arr))
        x_diag = x_base + N_arr * esc_n * tx
        y_diag = y_base + N_arr * esc_n * ty

        ax.fill(x_diag, y_diag, alpha=0.3, color='green')
        ax.plot(x_diag, y_diag, 'g-', linewidth=1.5)

        ax.text(xi + N_arr[0]*esc_n*tx*1.4, yi + N_arr[0]*esc_n*ty*1.4,
                f'{N_arr[0]:+.2f}', fontsize=8, ha='center', va='center', color='green', fontweight='bold')
        ax.text(xj + N_arr[-1]*esc_n*tx*1.4, yj + N_arr[-1]*esc_n*ty*1.4,
                f'{N_arr[-1]:+.2f}', fontsize=8, ha='center', va='center', color='green', fontweight='bold')

    for tag in [1, 2, 30, 4, 5]:
        ax.annotate(nombres_nodo[tag], (coords[tag][0], coords[tag][1]),
                    textcoords="offset points", xytext=(0, 12), ha='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.text(0.02, 0.02, 'Axial en tf (+ = tension, - = compresion)\nValores en extremos de cada elemento',
            transform=ax.transAxes, fontsize=9, va='bottom', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

    # ========================================
    # 5. DESPLAZAMIENTOS VERTICALES
    # ========================================
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_aspect('equal')
    ax.set_title('Desplazamientos Verticales Uy\n(Puerto a dos aguas - Rotula en C, q=3 tf/m)', fontsize=13, fontweight='bold')

    esc_d = 30

    for tag_i in range(len(tags_ord)-1):
        t1, t2 = tags_ord[tag_i], tags_ord[tag_i+1]
        if t1 == 30 and t2 == 31:
            continue
        x1, y1 = coords[t1]
        x2, y2 = coords[t2]
        ax.plot([x1, x2], [y1, y2], 'k--', linewidth=1, alpha=0.4)

    for tag in elements:
        ni, nj = elements[tag]
        xi, yi = coords[ni]
        xj, yj = coords[nj]
        n_pts = 50
        x_base = np.linspace(xi, xj, n_pts)
        y_base = np.linspace(yi, yj, n_pts)
        uy_interp = np.linspace(nodos[ni]['uy'], nodos[nj]['uy'], n_pts)

        x_def = x_base
        y_def = y_base + uy_interp * esc_d

        ax.plot(x_base, y_base, 'k-', linewidth=0.5, alpha=0.3)
        ax.plot(x_def, y_def, 'm-', linewidth=2, label='Deformada Uy' if tag==1 else '')
        ax.fill_between(x_def, y_base, y_def, alpha=0.2, color='magenta')

    for tag in [1, 2, 30, 4, 5]:
        uy_mm = nodos[tag]['uy'] * 1000
        ax.annotate(f'{nombres_nodo[tag]}\nUy={uy_mm:+.1f}mm',
                    (coords[tag][0], coords[tag][1]),
                    textcoords="offset points", xytext=(0, 15), ha='center', fontsize=9, fontweight='bold')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.text(0.02, 0.02, f'Deformada Uy amplificada x{esc_d}\n(Valores en mm)',
            transform=ax.transAxes, fontsize=9, va='bottom', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

    # ========================================
    # 6. TABLA DE RESULTADOS
    # ========================================
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.axis('off')
    ax.set_title('Resumen de Resultados\nPuerto a dos aguas - Rotula en C, q=3 tf/m proyectada', fontsize=14, fontweight='bold', pad=20)

    # Tabla de desplazamientos
    table_data = [['Nodo', 'Ux (mm)', 'Uy (mm)', 'Rot (rad)']]
    for tag in [1, 2, 30, 31, 4, 5]:
        n = nombres_nodo[tag]
        ux = nodos[tag]['ux']*1000
        uy = nodos[tag]['uy']*1000
        rz = nodos[tag]['rz']
        table_data.append([n, f'{ux:+.3f}', f'{uy:+.3f}', f'{rz:+.6f}'])

    table1 = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                      cellLoc='center', loc='upper left',
                      bbox=[0.0, 0.62, 0.55, 0.3])
    table1.auto_set_font_size(False)
    table1.set_fontsize(9)
    for key, cell in table1.get_celld().items():
        if key[0] == 0:
            cell.set_facecolor('#4472C4')
            cell.set_text_props(color='white', fontweight='bold')

    ax.text(0.0, 0.95, 'Desplazamientos Nodales', fontsize=11, fontweight='bold', transform=ax.transAxes)

    # Tabla de fuerzas
    ftable_data = [['Elem', 'Ext', 'Fx (tf)', 'Fy (tf)', 'Mz (tf.m)']]
    for tag in elements:
        ni, nj = elements[tag]
        f = fuerzas[tag]
        ftable_data.append([nombres_elem[tag], f'Nodo {ni}', f'{f[0]:+.3f}', f'{f[1]:+.3f}', f'{f[2]:+.3f}'])
        ftable_data.append(['', f'Nodo {nj}', f'{f[3]:+.3f}', f'{f[4]:+.3f}', f'{f[5]:+.3f}'])

    table2 = ax.table(cellText=ftable_data[1:], colLabels=ftable_data[0],
                      cellLoc='center', loc='upper right',
                      bbox=[0.58, 0.35, 0.42, 0.57])
    table2.auto_set_font_size(False)
    table2.set_fontsize(8)
    for key, cell in table2.get_celld().items():
        if key[0] == 0:
            cell.set_facecolor('#C0504D')
            cell.set_text_props(color='white', fontweight='bold')

    ax.text(0.58, 0.95, 'Fuerzas Internas (globales)', fontsize=11, fontweight='bold', transform=ax.transAxes)

    # Tabla de reacciones
    rtable_data = [['Apoyo', 'Rx (tf)', 'Ry (tf)', 'Mz (tf.m)']]
    for tag in [1, 5]:
        fx_s, fy_s, mz_s = 0, 0, 0
        for etag, (eni, enj) in elements.items():
            f = fuerzas[etag]
            if eni == tag:
                fx_s += f[0]; fy_s += f[1]; mz_s += f[2]
            elif enj == tag:
                fx_s += f[3]; fy_s += f[4]; mz_s += f[5]
        rtable_data.append([nombres_nodo[tag], f'{fx_s:+.3f}', f'{fy_s:+.3f}', f'{mz_s:+.3f}'])

    table3 = ax.table(cellText=rtable_data[1:], colLabels=rtable_data[0],
                      cellLoc='center', loc='lower left',
                      bbox=[0.0, 0.05, 0.55, 0.22])
    table3.auto_set_font_size(False)
    table3.set_fontsize(9)
    for key, cell in table3.get_celld().items():
        if key[0] == 0:
            cell.set_facecolor('#548235')
            cell.set_text_props(color='white', fontweight='bold')

    ax.text(0.0, 0.30, 'Reacciones', fontsize=11, fontweight='bold', transform=ax.transAxes)

    # Info
    ax.text(0.58, 0.25, f'Carga: q = 3.0 tf/m (proyectada vertical)\n'
                        f'Carga total: 39.0 tf\n'
                        f'SumFy = 39.000 tf\n'
                        f'SumFx = 0.000 tf\n'
                        f'Momento en rotula C = 0.000 tf.m',
            fontsize=10, transform=ax.transAxes, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close()

print("PDF generado: diagramas_puerto.pdf")
print("Ruta completa: C:\\Users\\edfev\\OneDrive\\Desktop\\UAndes\\noveno semestre\\Metodos computacionales en obras civiles\\Tarea 2\\diagramas_puerto.pdf")
