import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter

# ─── Materiais ────────────────────────────────────────────
materiais = {
    "Alumínio":  {"sigma": 3.5e7, "mu_r": 1.0,    "cor": "#3498db"},
    "Cobre":     {"sigma": 5.8e7, "mu_r": 1.0,    "cor": "#e67e22"},
    "Ferro":     {"sigma": 1.0e7, "mu_r": 1000.0, "cor": "#e74c3c"},
}

# ─── Parâmetros comuns ────────────────────────────────────
Lx, Ly = 0.15, 0.15
Nx, Ny = 100, 100
dx, dy = Lx/Nx, Ly/Ny
mu0    = 4 * np.pi * 1e-7
f      = 60.0
omega  = 2 * np.pi * f
B0     = 0.1
x      = np.linspace(0, Lx, Nx)
y      = np.linspace(0, Ly, Ny)
xg     = np.linspace(0, Lx*100, Nx)
yg     = np.linspace(0, Ly*100, Ny)

# ─── Simulação por material ───────────────────────────────
resultados = {}

for nome, props in materiais.items():
    sigma = props["sigma"]
    mu    = mu0 * props["mu_r"]
    alpha = 1.0 / (mu * sigma)
    dt    = 0.2 * mu * sigma * min(dx, dy)**2
    T     = 1.0 / f
    n_steps = int(T / dt)

    Bz      = np.zeros((Nx, Ny))
    losses  = np.zeros((Nx, Ny))
    Bz_hist = np.zeros((Nx, n_steps))
    snap_Jx = snap_Jy = None
    snap_ok = False
    t = 0.0

    for step in range(n_steps):
        phase = omega * t
        lap = (
            (np.roll(Bz,-1,0) + np.roll(Bz,1,0) - 2*Bz) / dx**2 +
            (np.roll(Bz,-1,1) + np.roll(Bz,1,1) - 2*Bz) / dy**2
        )
        Bz += dt * (alpha * lap + B0 * omega * np.cos(phase))
        Bext = B0 * np.sin(phase)
        Bz[0,:] = Bz[-1,:] = Bz[:,0] = Bz[:,-1] = Bext

        Jx = -np.gradient(Bz, dy, axis=1) / mu
        Jy =  np.gradient(Bz, dx, axis=0) / mu
        losses += (Jx**2 + Jy**2) * dt / sigma
        Bz_hist[:, step] = Bz[:, Ny//2]

        if not snap_ok and abs((phase % (2*np.pi)) - np.pi/2) < omega*dt*2:
            snap_Jx, snap_Jy = Jx.copy(), Jy.copy()
            snap_ok = True
        t += dt

    losses_avg = losses / T
    J_mag      = np.sqrt(snap_Jx**2 + snap_Jy**2) / 1e6
    Bz_perfil  = (Bz_hist.max(axis=1) - Bz_hist.min(axis=1)) / 2 * 1000
    delta      = np.sqrt(2 / (omega * mu * sigma))
    phi_lag    = (Lx/4) / delta
    lag_ms     = phi_lag / omega * 1000
    espessura  = 0.01
    P_total    = losses_avg.mean() * Lx * Ly * espessura

    resultados[nome] = {
        "delta":      delta,
        "lag_ms":     lag_ms,
        "phi_lag":    phi_lag,
        "J_mag":      J_mag,
        "J_max":      J_mag.max(),
        "J_mean":     J_mag.mean(),
        "losses_avg": losses_avg,
        "P_max":      losses_avg.max(),
        "P_mean":     losses_avg.mean(),
        "P_total":    P_total,
        "Bz_perfil":  Bz_perfil,
        "snap_Jx":    snap_Jx,
        "snap_Jy":    snap_Jy,
        "cor":        props["cor"],
        "sigma":      sigma,
        "mu_r":       props["mu_r"],
        "n_steps":    n_steps,
    }
    print(f"  ✓ {nome} simulado  ({n_steps} passos, dt={dt*1e6:.3f} μs)")

# ─── CLI ──────────────────────────────────────────────────
print("\n" + "=" * 62)
print("  COMPARAÇÃO DE MATERIAIS — CORRENTES DE FOUCAULT (60 Hz)")
print("=" * 62)
hdr = f"  {'Grandeza':<28} {'Alumínio':>10} {'Cobre':>10} {'Ferro':>10}"
print(hdr)
print("  " + "-" * 60)

linhas = [
    ("Condutividade [S/m]",  "sigma",   "{:.2e}"),
    ("Permeabilidade (μᵣ)",  "mu_r",    "{:.0f}"),
    ("Skin depth δ [mm]",    "delta",   "{:.2f}",   1000),
    ("Defasagem [ms]",       "lag_ms",  "{:.3f}"),
    ("Defasagem [°]",        "phi_lag", "{:.1f}",   180/np.pi),
    ("|J| máximo [kA/m²]",  "J_max",   "{:.2f}",   1e3),
    ("|J| médio [kA/m²]",   "J_mean",  "{:.2f}",   1e3),
    ("Razão pico/média",     None,      "{:.2f}"),
    ("P máxima [kW/m³]",    "P_max",   "{:.4f}",   1e-3),
    ("P média [kW/m³]",     "P_mean",  "{:.4f}",   1e-3),
    ("P total [mW]",         "P_total", "{:.4f}",   1e3),
]

for linha in linhas:
    label = linha[0]
    key   = linha[1]
    fmt   = linha[2]
    scale = linha[3] if len(linha) > 3 else 1.0

    vals = []
    for nome in materiais:
        r = resultados[nome]
        if key is None:
            v = r["J_max"] / r["J_mean"]
        else:
            v = r[key] * scale
        vals.append(fmt.format(v))

    print(f"  {label:<28} {vals[0]:>10} {vals[1]:>10} {vals[2]:>10}")

print("=" * 62)

# ─── Figura 4×3 ───────────────────────────────────────────
nomes = list(materiais.keys())
fig   = plt.figure(figsize=(16, 14))
fig.suptitle(
    f'Correntes de Foucault — Comparação de Materiais (f = {f:.0f} Hz, B₀ = {B0*1000:.0f} mT)',
    fontsize=14, fontweight='bold'
)
gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.48, wspace=0.35,
                       left=0.07, right=0.97, top=0.94, bottom=0.06)

for col, nome in enumerate(nomes):
    r   = resultados[nome]
    cor = r["cor"]

    # ── Linha 0: Densidade de corrente ────────────────────
    ax = fig.add_subplot(gs[0, col])
    im = ax.imshow(r["J_mag"].T, origin='lower', aspect='equal',
                   extent=[0, Lx*100, 0, Ly*100], cmap='inferno')
    cb = plt.colorbar(im, ax=ax, pad=0.02)
    cb.set_label('|J| [MA/m²]', fontsize=8)
    strm = ax.streamplot(xg, yg, r["snap_Jx"].T, r["snap_Jy"].T,
                         color='white', linewidth=0.5, density=1.2, arrowsize=0.7)
    strm.lines.set_alpha(0.45)
    ax.set_title(f'{nome}\n|J| máx = {r["J_max"]:.4f} MA/m²', fontsize=10, color=cor)
    ax.set_xlabel('x [cm]', fontsize=8)
    ax.set_ylabel('y [cm]', fontsize=8)

    # ── Linha 1: Perdas Joule ─────────────────────────────
    ax2 = fig.add_subplot(gs[1, col])
    im2 = ax2.imshow((r["losses_avg"]/1e3).T, origin='lower', aspect='equal',
                     extent=[0, Lx*100, 0, Ly*100], cmap='hot')
    cb2 = plt.colorbar(im2, ax=ax2, pad=0.02)
    cb2.set_label('P [kW/m³]', fontsize=8)
    ax2.set_title(f'P total = {r["P_total"]*1000:.3f} mW', fontsize=10)
    ax2.set_xlabel('x [cm]', fontsize=8)
    ax2.set_ylabel('y [cm]', fontsize=8)

# ── Linha 2: Perfil de atenuação Bz — todos juntos ────────
ax3 = fig.add_subplot(gs[2, :])
B_teo_al  = B0 * np.exp(-x / resultados["Alumínio"]["delta"]) * 1000
B_teo_cu  = B0 * np.exp(-x / resultados["Cobre"]["delta"])    * 1000
B_teo_fe  = B0 * np.exp(-x / resultados["Ferro"]["delta"])    * 1000

for nome, B_teo in zip(nomes, [B_teo_al, B_teo_cu, B_teo_fe]):
    r   = resultados[nome]
    cor = r["cor"]
    ax3.plot(x*100, B_teo, color=cor, lw=1.8, linestyle='--', alpha=0.7)
    ax3.plot(x*100, r["Bz_perfil"], color=cor, lw=2.2,
             label=f'{nome}  (δ={r["delta"]*1000:.1f} mm)')
    ax3.axvline(r["delta"]*100, color=cor, lw=1.0, linestyle=':', alpha=0.6)

ax3.axhline(B0*np.exp(-1)*1000, color='gray', lw=1.0, linestyle=':',
            label=f'B₀/e = {B0*np.exp(-1)*1000:.1f} mT')
ax3.set_title(r'Perfil de Atenuação de $B_z$ — Simulado (linha sólida) vs Teórico (tracejado)',
              fontsize=11)
ax3.set_xlabel('x [cm]', fontsize=10)
ax3.set_ylabel('|Bz| [mT]', fontsize=10)
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, Lx*100)
ax3.set_ylim(bottom=0)

# ── Linha 3: Evolução temporal — todos juntos ──────────────
ax4 = fig.add_subplot(gs[3, :])
t_arr = np.linspace(0, 2/f, 500)
B_ext_t = B0 * np.sin(omega * t_arr)
ax4.plot(t_arr*1000, B_ext_t*1000, color='gray', lw=2.0,
         linestyle='--', label='Campo externo $B_{ext}$', alpha=0.7)

for nome in nomes:
    r   = resultados[nome]
    cor = r["cor"]
    pen = Lx / 4
    att = np.exp(-pen / r["delta"])
    B_int = B0 * att * np.sin(omega * t_arr - pen / r["delta"])
    defasagem_graus = np.degrees(r["phi_lag"])
    ax4.plot(t_arr*1000, B_int*1000, color=cor, lw=2.2,
             label=f'{nome}  (Δφ={defasagem_graus:.1f}°, att={att:.3f})')
    ax4.fill_between(t_arr*1000, B_int*1000, alpha=0.08, color=cor)

ax4.set_title(r'Evolução Temporal — Defasagem e Atenuação por Material (x = 3.75 cm)',
              fontsize=11)
ax4.set_xlabel('t [ms]', fontsize=10)
ax4.set_ylabel('B [mT]', fontsize=10)
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, t_arr[-1]*1000)

plt.savefig('foucault_comparacao_materiais.png', dpi=180, bbox_inches='tight')
plt.show()
