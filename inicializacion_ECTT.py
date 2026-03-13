import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# PROYECTO ECTT - Gemelo Digital de Ensayo de Conductividad Térmica Transitoria
# Autor: [Tu Nombre]
# Descripción: Modelo 1D de transferencia de calor con cambio de fase y 
# convección natural (Rayleigh-Bénard en medios porosos).
# =============================================================================

# -----------------------------------------------------------------------------
# 1. PARÁMETROS GEOMÉTRICOS DEL ENSAYO
# -----------------------------------------------------------------------------
H = 0.040           # Altura de la muestra (m) -> 40 mm
D = 0.030           # Diámetro de la muestra (m) -> 30 mm
A_base = np.pi * (D / 2)**2  # Área de transferencia de calor (m^2)

# -----------------------------------------------------------------------------
# 2. PROPIEDADES DE LA ESPUMA METÁLICA (Aluminio)
# -----------------------------------------------------------------------------
PPI = 10            # Poros por pulgada
eps = 0.935         # Porosidad macroscópica (93.5%)
k_al = 200.0        # Conductividad térmica del aluminio (W/m·K)
rho_al = 2700.0     # Densidad del aluminio sólido (kg/m^3)
cp_al = 900.0       # Calor específico del aluminio (J/kg·K)

# Fórmulas geométricas (basadas en Ec. 156 y 157 de la bibliografía)
d_p = 0.0254 / PPI  # Diámetro del poro (m)
d_f = d_p * np.sqrt((1 - eps) / (3 * np.pi)) # Diámetro de la fibra/ligamento (m)

# Permeabilidad (Modelo de Carman-Kozeny adaptado)
K_perm = (d_p**2 * eps**3) / (150 * (1 - eps)**2) # Permeabilidad (m^2)

# -----------------------------------------------------------------------------
# 3. PROPIEDADES DEL MATERIAL DE CAMBIO DE FASE (PCM - Parafina tipo RT)
# -----------------------------------------------------------------------------
# Temperaturas de Fusión
T_solidus = 41.0    # Temperatura inicio fusión (ºC)
T_liquidus = 44.0   # Temperatura fin fusión (ºC)

# Propiedades Térmicas
L_f = 241200.0      # Calor latente de fusión (J/kg) -> Aprox. 67 Wh/kg
k_pcm_s = 0.2       # Conductividad térmica PCM sólido (W/m·K)
k_pcm_l = 0.2       # Conductividad térmica PCM líquido (W/m·K)
cp_pcm = 2000.0     # Calor específico PCM (J/kg·K)

# Propiedades Fluidodinámicas (Líquido)
rho_pcm_s = 880.0   # Densidad PCM sólido (kg/m^3)
rho_pcm_l = 760.0   # Densidad PCM líquido (kg/m^3)
mu_f = 0.003        # Viscosidad dinámica aproximada (Pa·s o kg/m·s)
nu_f = mu_f / rho_pcm_l  # Viscosidad cinemática (m^2/s)

# Coeficiente de expansión volumétrica (Aproximación de Boussinesq)
beta = (rho_pcm_s - rho_pcm_l) / (rho_pcm_l * (T_liquidus - T_solidus))

# Difusividad térmica del fluido (alfa_f)
alfa_f = k_pcm_l / (rho_pcm_l * cp_pcm) 
g = 9.81            # Gravedad (m/s^2)

# -----------------------------------------------------------------------------
# 4. CONDICIONES OPERATIVAS Y DE CONTORNO (El Hardware)
# -----------------------------------------------------------------------------
# Parámetros de la placa calefactora
V_fuente = 5.02     # Voltaje regulado en la fuente (V)
I_fuente = 0.42     # Intensidad fijada en modo CC (A)
Potencia_neta = V_fuente * I_fuente # Potencia real inyectada (~2.1 W)

# Flujo de calor constante (q'') en la base
q_flux = Potencia_neta / A_base     # W/m^2

# Condiciones Ambientales
T_amb = 20.0        # Temperatura ambiente del laboratorio (ºC)
U_perdidas = 0.5    # Coeficiente global de pérdidas estimadas (W/K)

# -----------------------------------------------------------------------------
# 5. DISCRETIZACIÓN ESPACIO-TEMPORAL (La Malla FDM)
# -----------------------------------------------------------------------------
n_nodos = 40        # Número de divisiones (1 milímetro por nodo)
dx = H / n_nodos    # Paso espacial (m)
dt = 0.5            # Paso temporal (segundos) -> Frecuencia de muestreo teórica
tiempo_sim = 3600   # Tiempo total de simulación (1 hora)
n_pasos = int(tiempo_sim / dt)

# Vector de Temperaturas Inicial (Todo empieza a temperatura ambiente)
T = np.ones(n_nodos) * T_amb
# Vector para guardar el historial del nodo del núcleo (ej. nodo central, índice 20)
historial_T_nucleo = []

print("Inicialización completada con éxito.")
print(f"Potencia neta inyectada: {Potencia_neta:.2f} W")
print(f"Permeabilidad calculada (K): {K_perm:.2e} m^2")
