import numpy as np
import matplotlib.pyplot as plt

Ms = np.array([1, 2, 4, 8])
h = 2.0 / Ms  # 高度方向特征尺寸

# 节点 3 Y 位移（绝对值）
disps = np.array([0.003930, 0.005254, 0.005991, 0.006206])

# Timoshenko 梁理论解（考虑剪切变形）
E, nu, L, H, t, P = 1000, 0, 2, 2, 1, 1
I = t * H**3 / 12  # 2/3
G = E / (2 * (1 + nu))  # 500
kappa = 5 / 6  # 矩形截面剪切修正系数
A = t * H  # 2

delta_bending = P * L**3 / (3 * E * I)  # 0.004
delta_shear = P * L / (kappa * G * A)  # 0.0024
exact_timoshenko = delta_bending + delta_shear  # 0.0064

error = np.abs(disps - exact_timoshenko) / exact_timoshenko

# 拟合收敛阶
p = np.polyfit(np.log(h), np.log(error), 1)[0]

plt.figure(figsize=(8, 6))
plt.loglog(h,
           error,
           'gs-',
           linewidth=2,
           markersize=10,
           markerfacecolor='lime',
           markeredgecolor='darkgreen',
           label=f'Height Refined T3 (p≈{p:.2f})')

# 参考线
h_ref = np.array([h[0], h[-1]])
plt.loglog(h_ref, (h_ref / h[0])**1 * error[0],
           'k--',
           linewidth=1.5,
           label='O(h)')
plt.loglog(h_ref, (h_ref / h[0])**2 * error[0],
           'k:',
           linewidth=1.5,
           label='O(h²)')

plt.xlabel('Element size in height $h_y = H/M$', fontsize=13)
plt.ylabel(
    'Relative error $|\\delta_h - \\delta_{Timoshenko}| / \\delta_{Timoshenko}$',
    fontsize=12)
plt.title('Convergence by Height Refinement ($L=2, H=2$, N=2 fixed)\n' +
          r'T3 Element vs. Timoshenko Beam Solution $\delta_{exact}=0.0064$',
          fontsize=11)
plt.legend(fontsize=11, loc='upper left')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('convergence_height_timoshenko.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Timoshenko exact: {exact_timoshenko:.4f}")
print(f"Observed convergence rate: {p:.2f}")
for i, M in enumerate(Ms):
    print(f"M={M}: disp={disps[i]:.5f}, err={error[i]*100:.1f}%")
