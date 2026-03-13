import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# GEMELO DIGITAL ECTT - MODELO UNIFICADO (LISTO PARA EJECUTAR)
# =============================================================================

print("Inicializando variables y parámetros geométricos...")

# 1. PARÁMETROS GEOMÉTRICOS
H = 0.040           # Altura de la muestra (40 mm)
D = 0.030           # Diámetro de la muestra (30 mm)
A_base = np.pi * (D / 2)**2

# 2. PROPIEDADES DE LA ESPUMA (Aluminio)
PPI = 10
eps = 0.935         # Porosidad
k_al = 200.0
rho_al = 2700.0
cp_al = 900.0

d_p = 0.0254 / PPI
K_perm = (d_p**2 * eps**3) / (150 * (1 - eps)**2) # Permeabilidad

# 3. PROPIEDADES DEL PCM
T_solidus = 41.0
T_liquidus = 44.0
L_f = 241200.0
k_pcm_s = 0.2
k_pcm_l = 0.2
cp_pcm = 2000.0
rho_pcm_s = 880.0
rho_pcm_l = 760.0

mu_f = 0.003
nu_f = mu_f / rho_pcm_l
beta = (rho_pcm_s - rho_pcm_l) / (rho_pcm_l * (T_liquidus - T_solidus))
alfa_f = k_pcm_l / (rho_pcm_l * cp_pcm)
g = 9.81

# 4. CONDICIONES OPERATIVAS
V_fuente = 5.02
I_fuente = 0.42
Potencia_neta = V_fuente * I_fuente
q_flux = Potencia_neta / A_base
T_amb = 20.0

# 5. DISCRETIZACIÓN (MALLA Y TIEMPO)
n_nodos = 40
dx = H / n_nodos
dt = 0.1            # ¡Reducido a 0.1s para garantizar estabilidad matemática absoluta!
tiempo_sim = 3600   # 1 hora
n_pasos = int(tiempo_sim / dt)

T = np.ones(n_nodos) * T_amb

# 6. PROPIEDADES EFECTIVAS ESTÁTICAS
k_eff_solido = (1 - eps) * k_al + eps * k_pcm_s
k_eff_liquido_puro = (1 - eps) * k_al + eps * k_pcm_l
rho_cp_solido = (1 - eps) * rho_al * cp_al + eps * rho_pcm_s * cp_pcm
rho_cp_liquido = (1 - eps) * rho_al * cp_al + eps * rho_pcm_l * cp_pcm

T_old = np.copy(T)
T_new = np.zeros(n_nodos)
alfa = np.zeros(n_nodos)
k_eff_dinamico = np.ones(n_nodos) * k_eff_solido
Cp_app = np.zeros(n_nodos)

tiempo_plot = []
T_base_plot = []
T_nucleo_plot = []

print("Simulando 1 hora de ensayo. Por favor, espera unos segundos...")

# =============================================================================
# BUCLE TRANSIORIO PRINCIPAL
# =============================================================================
for step in range(n_pasos):
    t_actual = step * dt
    
    # 1. ESTADO FÍSICO (Cambio de Fase)
    for i in range(n_nodos):
        if T_old[i] < T_solidus:
            alfa[i] = 0.0
            Cp_app[i] = rho_cp_solido
        elif T_old[i] > T_liquidus:
            alfa[i] = 1.0
            Cp_app[i] = rho_cp_liquido
        else:
            alfa[i] = (T_old[i] - T_solidus) / (T_liquidus - T_solidus)
            calor_latente_vol = (eps * rho_pcm_s * L_f) / (T_liquidus - T_solidus)
            Cp_app[i] = rho_cp_solido + calor_latente_vol

    # 2. CONVECCIÓN (Rayleigh)
    H_liquido = np.sum(alfa) * dx
    if T_old[0] > T_liquidus and H_liquido > 0:
        Ra_K = (g * beta * (T_old[0] - T_liquidus) * K_perm * H_liquido) / (nu_f * alfa_f)
        Nu = (Ra_K / 40) ** 0.5 if Ra_K >= 40 else 1.0
    else:
        Nu = 1.0

    for i in range(n_nodos):
        if alfa[i] == 1.0:
            k_eff_dinamico[i] = k_eff_liquido_puro * Nu
        else:
            k_eff_dinamico[i] = k_eff_solido

    # 3. TRANSFERENCIA DE CALOR (Ley de Fourier)
    # Base (Nodo 0)
    T_new[0] = T_old[0] + (dt / Cp_app[0]) * ((q_flux / dx) - k_eff_dinamico[0] * (T_old[0] - T_old[1]) / dx**2)
    
    # Medio (Nodos 1 a 38)
    for i in range(1, n_nodos - 1):
        q_in = k_eff_dinamico[i] * (T_old[i-1] - T_old[i]) / dx**2
        q_out = k_eff_dinamico[i] * (T_old[i] - T_old[i+1]) / dx**2
        T_new[i] = T_old[i] + (dt / Cp_app[i]) * (q_in - q_out)
        
    # Techo (Nodo 39 - asilado)
    T_new[n_nodos-1] = T_old[n_nodos-1] + (dt / Cp_app[n_nodos-1]) * (k_eff_dinamico[n_nodos-1] * (T_old[n_nodos-2] - T_old[n_nodos-1]) / dx**2)

    T_old = np.copy(T_new)
    
    # Guardamos datos para la gráfica cada 10 segundos
    if step % int(10 / dt) == 0:
        tiempo_plot.append(t_actual / 60)
        T_base_plot.append(T_new[0])
        T_nucleo_plot.append(T_new[20]) # Nodo a 20mm de altura

print("¡Cálculos terminados! Generando gráfica...")

# =============================================================================
# GRÁFICA DE RESULTADOS
# =============================================================================
plt.figure(figsize=(10, 6))
plt.plot(tiempo_plot, T_base_plot, label='T Base (0 mm)', color='red', linewidth=2)
plt.plot(tiempo_plot, T_nucleo_plot, label='T Núcleo (20 mm)', color='blue', linewidth=2)

plt.axhline(y=T_solidus, color='gray', linestyle='--', label='Inicio Fusión (41ºC)')
plt.axhline(y=T_liquidus, color='black', linestyle='--', label='Fin Fusión (44ºC)')

plt.title('Gemelo Digital ECTT: Curvas de Calentamiento')
plt.xlabel('Tiempo (minutos)')
plt.ylabel('Temperatura (ºC)')
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()
