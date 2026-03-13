import numpy as np
import matplotlib.pyplot as plt

# --- PRE-CÁLCULOS ESTÁTICOS ---
# Conductividad efectiva estática (Muestra totalmente sólida)
# Usamos un modelo simplificado de mezcla ponderada para la espuma y el PCM
k_eff_solido = (1 - eps) * k_al + eps * k_pcm_s  # Aprox. 13 W/mK
k_eff_liquido_puro = (1 - eps) * k_al + eps * k_pcm_l

# Capacidad calorífica volumétrica efectiva (rho * cp)
rho_cp_solido = (1 - eps) * rho_al * cp_al + eps * rho_pcm_s * cp_pcm
rho_cp_liquido = (1 - eps) * rho_al * cp_al + eps * rho_pcm_l * cp_pcm

# Preparar arrays para el bucle
T_old = np.copy(T)        # Temperaturas en el instante actual
T_new = np.zeros(n_nodos) # Temperaturas en el instante siguiente
alfa = np.zeros(n_nodos)  # Fracción líquida en cada nodo
k_eff_dinamico = np.ones(n_nodos) * k_eff_solido 
Cp_app = np.zeros(n_nodos)

# Variables para graficar
tiempo_plot = []
T_base_plot = []
T_nucleo_plot = []  # Nodo central (mitad de la probeta)

print("Iniciando la simulación termodinámica... (Esto puede tardar unos segundos)")

# =============================================================================
# BUCLE DE TIEMPO PRINCIPAL
# =============================================================================
for step in range(n_pasos):
    t_actual = step * dt
    
    # -------------------------------------------------------------------------
    # FASE 1: EVALUAR EL ESTADO FÍSICO DE CADA NODO (Cambio de fase)
    # -------------------------------------------------------------------------
    for i in range(n_nodos):
        if T_old[i] < T_solidus:
            # Estado Sólido
            alfa[i] = 0.0
            Cp_app[i] = rho_cp_solido
        elif T_old[i] > T_liquidus:
            # Estado Líquido
            alfa[i] = 1.0
            Cp_app[i] = rho_cp_liquido
        else:
            # Zona pastosa (Transición de fase)
            alfa[i] = (T_old[i] - T_solidus) / (T_liquidus - T_solidus)
            # Método de Capacidad Calorífica Aparente (absorbe el calor latente)
            calor_latente_volumetrico = (eps * rho_pcm_s * L_f) / (T_liquidus - T_solidus)
            Cp_app[i] = rho_cp_solido + calor_latente_volumetrico

    # -------------------------------------------------------------------------
    # FASE 2: CÁLCULO DE LA CONVECCIÓN (Números de Rayleigh y Nusselt)
    # -------------------------------------------------------------------------
    # ¿Cuánto líquido hay en total? (La "pista de despegue" H_liquido)
    H_liquido = np.sum(alfa) * dx
    
    # Solo calculamos convección si la base está más caliente que la fusión
    if T_old[0] > T_liquidus and H_liquido > 0:
        Delta_T_conv = T_old[0] - T_liquidus
        
        # Número de Rayleigh de Darcy
        Ra_K = (g * beta * Delta_T_conv * K_perm * H_liquido) / (nu_f * alfa_f)
        
        # Criterio de Horton-Rogers-Lapwood (Ra_c = 40)
        if Ra_K >= 40:
            Nu = (Ra_K / 40) ** 0.5  # Correlación clásica para medios porosos
        else:
            Nu = 1.0
    else:
        Ra_K = 0
        Nu = 1.0

    # Actualizar la conductividad de cada nodo según su estado
    for i in range(n_nodos):
        if alfa[i] == 1.0:
            # Si es líquido total, aplicamos la mejora convectiva de tu modelo
            k_eff_dinamico[i] = k_eff_liquido_puro * Nu
        else:
            k_eff_dinamico[i] = k_eff_solido

    # -------------------------------------------------------------------------
    # FASE 3: LEY DE FOURIER (Ecuación de Diferencias Finitas 1D)
    # -------------------------------------------------------------------------
    # Nodo 0: Condición de Contorno de Flujo de Calor Constante (La placa de 2.1W)
    T_new[0] = T_old[0] + (dt / Cp_app[0]) * ((q_flux / dx) - k_eff_dinamico[0] * (T_old[0] - T_old[1]) / dx**2)
    
    # Nodos intermedios: Conducción de calor hacia arriba
    for i in range(1, n_nodos - 1):
        q_in = k_eff_dinamico[i] * (T_old[i-1] - T_old[i]) / dx**2
        q_out = k_eff_dinamico[i] * (T_old[i] - T_old[i+1]) / dx**2
        T_new[i] = T_old[i] + (dt / Cp_app[i]) * (q_in - q_out)
        
    # Nodo final (Arriba): Condición de contorno adiabática o con pérdidas
    # (Asumimos bien aislado arriba para simplificar este código base)
    T_new[n_nodos-1] = T_old[n_nodos-1] + (dt / Cp_app[n_nodos-1]) * (k_eff_dinamico[n_nodos-1] * (T_old[n_nodos-2] - T_old[n_nodos-1]) / dx**2)

    # Actualizar temperaturas para el siguiente segundo
    T_old = np.copy(T_new)
    
    # Guardar datos cada 10 segundos para la gráfica (optimizar memoria)
    if step % (10 / dt) == 0:
        tiempo_plot.append(t_actual / 60) # Guardamos en minutos
        T_base_plot.append(T_new[0])
        T_nucleo_plot.append(T_new[20])   # Nodo 20 es el centro exacto (20mm)

print("¡Simulación terminada!")

# =============================================================================
# VISUALIZACIÓN DE RESULTADOS (Validación Experimental)
# =============================================================================
plt.figure(figsize=(10, 6))
plt.plot(tiempo_plot, T_base_plot, label='T_base (Placa Calefactora)', color='red', linewidth=2)
plt.plot(tiempo_plot, T_nucleo_plot, label='T_núcleo (Termopar a 20mm)', color='blue', linewidth=2)

# Marcar las isotermas de fusión
plt.axhline(y=T_solidus, color='gray', linestyle='--', alpha=0.7, label='Inicio Fusión (Solidus)')
plt.axhline(y=T_liquidus, color='black', linestyle='--', alpha=0.7, label='Fin Fusión (Liquidus)')

plt.title('Gemelo Digital ECTT: Evolución Térmica del Ensayo')
plt.xlabel('Tiempo (minutos)')
plt.ylabel('Temperatura (ºC)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()
