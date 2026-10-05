# Cooling-Tower

# 🌬️ Water Cooling / Air Dehumidification Tower Simulator

A Streamlit-based engineering application for modeling and analysis of Water Cooling Towers and Air Dehumidification Towers.

The tool calculates the required tower height and generates graphical profiles for heat and mass transfer operations using operating line and saturation curve concepts.

---

## 🚀 Features

✅ Cooling Tower / Dehumidification Tower Design

✅ Tower Height Calculation

✅ Saturation Curve Visualization

✅ Operating Line Analysis

✅ Tie-Line Representation

✅ Humidity Profile Along Tower Height

✅ Liquid Temperature Profile

✅ Gas Temperature Profile

✅ Interactive Parameter Inputs

✅ Real-Time Streamlit Dashboard

---

## 📷 Application Screenshot

### Main Dashboard

screenshots/1.png

---

## ⚙️ Input Parameters

The model accepts the following inputs:

| Parameter | Description |
| ---------- | ---------- |
| TL2 (°C) | Outlet Liquid Temperature |
| TL1 (°C) | Inlet Liquid Temperature |
| TG1 (°C) | Inlet Gas Temperature |
| H1 (kg/kg) | Inlet Humidity Ratio |
| L/S | Liquid-to-Gas Flow Ratio |
| G'/S | Gas Flow Parameter |
| kYa | Volumetric Mass Transfer Coefficient |
| hLa | Liquid Side Heat Transfer Coefficient |

---

## 📊 Outputs

The application computes:

- Tower Height (m)
- Humidity Distribution
- Gas Temperature Distribution
- Liquid Temperature Distribution
- Operating Line
- Saturation Curve
- Tie Lines

---

## 📈 Generated Graphs

### 1. Operation Line

Displays:

- Saturation Curve
- Operating Line
- Tie Lines

Useful for visualizing heat and mass transfer driving forces.

### 2. Humidity Profile

Shows variation of humidity ratio through the tower height.

### 3. Liquid Temperature Profile

Visualizes cooling of liquid along the tower.

### 4. Gas Temperature Profile

Tracks the change in gas temperature through tower height.

---

## 🛠️ Technology Stack

- Python 3.10+
- Streamlit
- NumPy
- SciPy
- Matplotlib
- Pandas

---

## 📦 Installation

Clone the repository:

```bash
git clone https://github.com/bhaskar08092000/Water-Cooling-Air-Dehumidification-Tower.git
```

Navigate to project directory:

```bash
cd Water-Cooling-Air-Dehumidification-Tower
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## ▶️ Run Application

```bash
streamlit run app.py
```

Or

```bash
python -m streamlit run app.py
```

---

## 🎯 Applications

- Chemical Engineering
- Thermal Engineering
- HVAC Systems
- Cooling Tower Design
- Air Dehumidification Systems
- Heat & Mass Transfer Studies
- Academic Research
- Process Plant Design

---

## 📚 Engineering Concepts Used

The application implements principles of:

- Heat Transfer
- Mass Transfer
- Psychrometrics
- Cooling Tower Theory
- Merkel-Type Analysis
- Differential Tower Calculations
- Gas-Liquid Contact Operations

---

## 👨‍💻 Author

**Bhaskar Kewalramani**

Senior Engineer

GitHub: https://github.com/bhaskar08092000

---

## License

This project is released under the MIT License.
