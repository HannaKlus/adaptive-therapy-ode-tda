
#-LIBRARIES-
import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams
from scipy.integrate import solve_ivp

def model(t, y, drug_on):
    t_plus, t_P,t_minus, psa = y
    if drug_on:
        K_plus = 0.5 * t_P
        K_P = 100
    else:
        K_plus = 1.5 * t_P
        K_P = 10000
    K_minus = 10000
    r_plus = 0.0027726
    r_P = 0.0034657
    r_minus = 0.0066542

    K_plus = max(K_plus, 1e-9)

    dt_plus =  r_plus*t_plus*(1-(1*t_plus + 0.7*t_P+0.8*t_minus)/K_plus)
    dt_P = r_P*t_P*(1-(0.4*t_plus+1*t_P+0.5*t_minus)/K_P)
    dt_minus = r_minus*t_minus*(1 - (0.6*t_plus + 0.9*t_P + 1*t_minus)/K_minus)



    sigmaPSA = 0.5#rozpad PSA
    dpsa = (t_plus+t_P+t_minus)-sigmaPSA*psa
    return [dt_plus, dt_P, dt_minus, dpsa]

ess_vals = np.array([6060.60606060606,7575.75757575758,1.e-09,27272.7272727273])
ess_psa = ess_vals[3]
PSA_zenith = ess_psa * 0.8
PSA_nadir = PSA_zenith * 0.4

def event_psa_drop(t, y, drug_on):
    return y[3] - PSA_nadir

event_psa_drop.terminal = True
event_psa_drop.direction = -1

def event_psa_recover(t, y, drug_on):
    return y[3] - PSA_zenith

event_psa_recover.terminal = True
event_psa_recover.direction = 1

initial_vals = ess_vals * 0.4
time_current = 0
time_max = 17500
drug_on = False
time_points = []
psa_data = []

while time_current < time_max:
    current_event = event_psa_drop if drug_on else event_psa_recover
    result = solve_ivp(model, [time_current, time_max], initial_vals, args = (drug_on,), events = current_event, max_step = 1, method = "BDF")

    time_points.extend(result.t)
    psa_data.extend(result.y[3])

    time_current = result.t[-1]
    initial_vals = result.y[:, -1]

    if result.status == 0:
        break
    drug_on = not drug_on


#TDA - Taken's theorem and Time delayed embedding
tau = 350
X_tda = []

for i in range(len(psa_data) - tau):
    x_now = psa_data[i]
    x_delayed = psa_data[i + tau]
    X_tda.append([x_now, x_delayed])
X_tda = np.array(X_tda)

window_size = 1200
step_size = 50
h1_persistence = []
time_axis = []

for start in range(0, len(psa_data) - tau - window_size, step_size):
    end = start + window_size
    window = X_tda[start:end]

    window = window[::4]

    result = ripser(window, maxdim=1)
    diagrams = result['dgms']

    if len(diagrams[1]) == 0:
        h1_persistence.append(0.0)

    else:
        temorary_persistence = []
        for pt in diagrams[1]:
            persistence = pt[1] - pt[0]
            temorary_persistence.append(persistence)
        h1_persistence.append(np.max(temorary_persistence))
    time_axis.append(time_points[start + window_size//2])


import matplotlib.pyplot as plt


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), gridspec_kw={'wspace': 0.2})

ax1.plot(time_axis, h1_persistence, color='#d62728', linewidth=2)
ax1.set_title('Topological Mutation Detector ($H_1$ Persistence)', fontsize=16, fontweight='bold', pad=15)
ax1.set_xlabel('Time [days]', fontsize=14)
ax1.set_ylabel('Max $H_1$ Persistence', fontsize=14)


ax2.plot(time_points, psa_data, color='#9467bd', linewidth=2, label='PSA Level')
ax2.set_title('PSA Tumor Dynamics', fontsize=16, fontweight='bold', pad=15)
ax2.set_xlabel('Time [days]', fontsize=14)
ax2.set_ylabel('PSA Concentration', fontsize=14)
ax2.legend(fontsize=12, loc='upper left')


for ax in [ax1, ax2]:
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.savefig('Detect_Patient_2.png')
plt.show()





