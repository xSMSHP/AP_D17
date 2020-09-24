import re
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings("ignore")

class Passiv:
    def __init__(self, val, name, base):
        self.U2 = []
        self.U1eff = []
        self.Ueff = []
        self.f_base = []
        self.f = []
        self.name = name
        
        for line in val:
            l = re.split(',', line)
            self.U2.append((float(l[0]) + float(l[1])) / 2 * 10 ** (-3))
            self.U1eff.append(56 * 10 ** (-3) * (float(l[2]) + float(l[3])) / 2)
            self.f.append(float(l[4]))
        
        for line in base:
            l = re.split(',', line)
            self.Ueff.append(56 * float(l[0]) * 10 ** (-3))
            self.f_base.append(float(l[1]))

    def plot(self):
        Umax = 0
        fmax = 0
        for i, U in enumerate(self.U2):
            fmax = self.f[i] if U > Umax and i > 0 else fmax
            Umax = U if U > Umax and i > 0 else Umax
            
        plt.figure()
        plt.vlines(fmax, 0, 5, colors = 'k')
        plt.plot(self.f, self.U2, 'yo', self.f, self.U1eff, 'bd')
        plt.plot(self.f_base, self.Ueff, 'r*')
        plt.legend([r'$U_{2s}(\nu)$',r'$U_{1,eff}$',r'$U_{eff}$ ohne passiven Schwingkreis'], bbox_to_anchor = (0.75, -0.15))
        plt.grid(True)
        plt.ylabel('Spannung U [V]')
        plt.xlabel(r'Frequenz $\nu$ [MHz]')
        plt.autoscale(True)
        plt.show()
        plt.savefig(self.name + '.png', dpi = 900)

class Probe:
    def __init__(self, file, name, f, d, cutoff, cutin, ma, mi):
        self.name = name
        self.d = d
        self.t = []
        self.U1 = []
        self.U2 = []
        self.dB = 0
        self.B = 0
        self.f = f
        self.c = cutoff
        self.ci = cutin
        self.ma = ma
        self.mi = mi
        
        txt = file
        for i, line in enumerate(txt):
            if i > 4:
                tx_line = line.replace('\tNAN', '').replace(',', '.').replace('\t', ',')
                l = re.split(',', tx_line)
                self.t.append(float(l[0]))
                self.U1.append(float(l[1]))
                self.U2.append(float(l[2]))
    
    def plot(self):
        coef1 = np.polyfit(self.U2[self.ci:int(len(self.U2) / 2) + self.d], self.U1[self.ci:int(len(self.U1) / 2) + self.d], 55)
        coef2 = np.polyfit(self.U2[int(len(self.U2) / 2) + self.d:len(self.U2) - self.c], self.U1[int(len(self.U1) / 2) + self.d:len(self.U1) - self.c], 55)
        poly1d_fn1 = np.poly1d(coef1) 
        poly1d_fn2 = np.poly1d(coef2) 
        
        mini = [100, 100]
        index = [0, 0]
        ground = [0, 0]
        g_count = [0, 0]
        
        for i, U in enumerate(self.U1[self.ci:int(len(self.U1) / 2) + self.d]):
            if U < mini[0]:
                mini[0] = U
                index[0] = i + self.ci
            if U > 0:
                ground[0] += U
                g_count[0] += 1
        for i, U in enumerate(self.U1[int(len(self.U1) / 2) + self.d:len(self.U1) - self.c]):
            if U < mini[1]:
                mini[1] = U
                index[1] = i + int(len(self.U1) / 2) + self.d
            if U > 0:
                ground[1] += U
                g_count[1] += 1

        pos = [(mini[0] - ground[0] / g_count[0]) / 2, (mini[1] - ground[1] / g_count[1]) / 2]
        x1 = (np.poly1d(coef1) - pos[0] - ground[0] / g_count[0]).roots
        x2 = (np.poly1d(coef2) - pos[1] - ground[0] / g_count[0]).roots
        out = []
        
        for x in x1:
            x = np.real_if_close(x)
            if x < self.ma and x > self.mi:
                if np.isreal(x):
                    out.append(x) 
        for x in x2:
            x = np.real_if_close(x)
            if x < self.ma and x > self.mi:
                if np.isreal(x):
                    out.append(x)
                    
        plt.figure()
        self.B = 0.674 * (self.U2[index[0]] + self.U2[index[1]]) / 2
        if len(out) == 4:       
            self.dB = 0.674 * (out[0] - out[1] + out[2] - out[3]) / 2
            plt.hlines(pos[0] + ground[0] / g_count[0], out[1], out[0], colors = 'purple')
            plt.hlines(pos[1] + ground[1] / g_count[1], out[3], out[2], colors = 'orange')
        else: 
            print(out)
        
        plt.plot(self.U2[self.ci:int(len(self.U2) / 2) + self.d], self.U1[self.ci:int(len(self.U1) / 2) + self.d], color = 'b', marker = ',', alpha = 0.6)
        plt.plot(self.U2[int(len(self.U2) / 2) + self.d:len(self.U2) - self.c], self.U1[int(len(self.U1) / 2) + self.d:len(self.U1) - self.c], color = 'r', marker = ',', alpha = 0.6)
        plt.plot(np.linspace(self.mi, self.ma, 100), poly1d_fn1(np.linspace(self.mi, self.ma, 100)), 'b-')
        plt.plot(np.linspace(self.mi, self.ma, 100), poly1d_fn2(np.linspace(self.mi, self.ma, 100)), 'r-')
        plt.xlabel(r'$U_A$ [V]')
        plt.ylabel(r'$U_B$ [V]')
        plt.legend([r'$U_A$ Messwerte', r'$U_B$ Messwerte', r'$U_A$ Polynom 55. Grades', r'$U_B$ Polynom 55. Grades'])
        plt.autoscale(True)
        plt.grid(True)
        plt.show()
        plt.savefig(self.name + '.png', dpi = 900)
    
a1_txt = open("auf1.txt", "r")
a1_num = 1
a1_base = []
a1_current = []
last_mark = 0

for line in a1_txt:
    last_mark = -2 if line == '-2\n' else -1 if line == '-1\n' else last_mark
    if last_mark == 0:
        a1_base.append(line)
    elif last_mark == -2:
        if line == '-2\n':
            continue
        a1_current.append(line)
    else:
        p = Passiv(a1_current, 'kond' + str(a1_num), a1_base)
        p.plot()
        a1_num += 1
        a1_current = []
        last_mark = -2

a2_txt = [open("auf2_31Hz.txt", "r"), open("auf2_35Hz.txt", "r"), open("auf2_40Hz.txt", "r"), open("auf2_45Hz.txt", "r"), open("auf2_50Hz.txt", "r"), open("auf2_55Hz.txt", "r"), open("auf2_60_9Hz.txt", "r"), open("auf2_65Hz.txt", "r"), open("auf2_70Hz.txt", "r"), open("auf2_73Hz.txt", "r")]
freq = [31, 35, 40, 45, 50, 55, 60.9, 65, 70, 73]
delta = [0, 0, -100, 0, 700, 0, 0, 0, 0, 0]
c = [0, 500, 0, 0, 0, 0, 0, 1100, 0, 0]
ci = [0, 0, 0, 0, 400, 300, 0, 0, 0, 400]
ma = [2.8, 3, 3.2, 3.2, 4, 4, 4.5, 4.3, 5, 5]
mi = [1, 1, 1.3, 1.5, 2, 2, 2.3, 2.5, 3, 3]
B0 = []
dB0 = []

for i, a2 in enumerate(a2_txt):
    pr = Probe(a2, str(int(freq[i])) + 'MHz', freq[i], delta[i], c[i], ci[i], ma[i], mi[i])
    pr.plot()
    B0.append(pr.B)
    dB0.append(pr.dB)

coefB0 = np.polyfit(B0, freq, 1)
poly1d_fnB0 = np.poly1d(coefB0) 
h = 6.626 * 10 ** (-34)
muhB = 9.274 * 10 ** (-24)

dist = 0

for i, B in enumerate(B0):
    dist += abs(freq[i] - poly1d_fnB0(B))    

print('Dist = ' + str(dist / len(B0)))
print(poly1d_fnB0)
print('g = ' + str(h / muhB * coefB0[0] * 10 ** 9))

plt.figure()
plt.plot(B0, freq, 'bo')
plt.plot(np.linspace(0,3), poly1d_fnB0(np.linspace(0,3)), 'b-')
plt.legend(['Aus Messwerten', 'Ausgleichsgerade'])
plt.xlabel(r'$B_0$ [mT]')
plt.ylabel(r'$\nu_0$ [MHz]')
plt.grid(True)
plt.autoscale(True)
plt.show()
plt.savefig('B0.png', dpi = 900)

dB0av = 0

for dB in dB0:
    dB0av += dB

print('dB0 = ' + str(dB0av / len(dB0)))

plt.figure()
plt.plot(freq, dB0, 'bo')
plt.xlabel(r'$\delta B_0$ [mT]')
plt.ylabel(r'$\nu_0$ [MHz]')
plt.grid(True)
plt.autoscale(True)
plt.show()
plt.savefig('dB0.png', dpi = 900)
        