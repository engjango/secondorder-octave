# Second Order System Parameters

This script randomly chooses two Script_2Ord * .m files present
in the same directory and calculates, for each one, the following
second order system parameters:

- Natural frequency (wn);
- Damping factor (zeta);
- Transfer Function [G (s)];
- Peak instant (Tp);
- Percentage surpassing (UP);
- Accommodation time (Ts);
- Rising time (Tr);

![Screenshot main](https://raw.githubusercontent.com/jndgomes/secondorder-octave/main/screenshot.png)

---

# Parâmetros de Sistemas de Segunda Ordem (Sistemas de Controle)

Este script escolhe aleatoriamente dois arquivos Script_2Ord*.m presentes
no mesmo diretório e realiza o cálculo, para cada um, dos seguintes 
parâmetros dos sistemas de segunda ordem, respectivamente:

- Frequência natural (wn); 
- Fator de amortecimento (zeta); 
- Função de Transferência [G(s)];
- Instante de pico (Tp);
- Ultrapassagem de percentual (UP);
- Tempo de Acomodação (Ts);
- Tempo de subida (Tr);

---

# Output

Some results obtained for files Sist_2Ord3.mat and Sist_2Ord5.mat.

```
>> script_sun
-----------------------------------------------------------
Arquivo 1: Sist_2Ord3.mat

Função de Transferência:
G(s)_1 = 0.53/(s^2 + 0.84s + 1.00).

Parâmetros:
Tp = 3.46 s
Ts = 8.38 s
Tr = 1.98s - 0.48s = 1.50 s
Wn =  1.00 rad/s
Zeta =  0.42
UP% = 23.38%

Parâmetros Re-Calculados:
Tp = 3.46 s
Ts = 9.52 s
Tr = 1.46 s
Wn =  1.00 rad/s
Zeta =  0.42
UP% = 23.38%
-----------------------------------------------------------

-----------------------------------------------------------
Arquivo 2: Sist_2Ord5.mat

Função de Transferência:
G(s)_2 = 1.31/(s^2 + 0.97s + 1.21).

Parâmetros:
Tp = 3.18 s
Ts = 7.59 s
Tr = 1.83s - 0.44s = 1.39 s
Wn =  1.10 rad/s
Zeta =  0.44
UP% = 21.46%

Parâmetros Re-Calculados:
Tp = 3.18 s
Ts = 8.27 s
Tr = 1.33 s
Wn =  1.10 rad/s
Zeta =  0.44
UP% = 21.46%
-----------------------------------------------------------
```

# Install Dependencies

## control package

If you haven't yet installed the control package, then type the following command in the GNU/Octave terminal window and press Enter.

```
>> pkg install -forge control
```

# Compile and Run

Once you have the dependencie (see above), open the file script_sun.m on GNU/Octave 
and go to `Execute > Save and Execute` (or just press F5).

# Dependencies

* control package.
