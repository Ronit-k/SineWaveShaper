# Two-Threshold Waveshaper — Triangle-to-Sine Conversion  
### EE380 – Digital Circuits Lab (IIT Kanpur)

This MATLAB script simulates and analyzes the conversion of a **±4 V triangular wave** into an **approximate sine wave** using a **two-threshold nonlinear waveshaper**.  

It demonstrates:

- Manual piecewise VTC design (as in IITK lab manual)
- Search-based optimal threshold selection
- THD calculation up to 10th harmonic (lab-style)
- Waveform & FFT visualization for input and output
- Option to enforce **symmetry constraint** `y = x` at ±t₁ (diode-like transfer curve)

---

## 📂 File

| File | Description |
|------|------------|
| `SineOptimizerLab6.m` | Main script for simulation, optimization & plotting |

---

## 🛠️ How to Run

Open MATLAB and execute:

```matlab
>> SineOptimizerLab6
