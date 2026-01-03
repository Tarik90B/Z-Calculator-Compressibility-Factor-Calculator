# Z Calculator – Compressibility Factor Calculator

**Z Calculator** is a Python/Tkinter application that calculates the **compressibility factor (Z)** for natural gas using the **AGA8 equation of state** (ISO 12213).  

This software is intended for educational and engineering reference purposes.

---

## Disclaimer

This software implements the AGA8 equation of state based on ISO 12213 (Natural gas — Calculation of compression factor).  

This is an independent software implementation. **ISO and the American Gas Association (AGA) do not endorse, certify, or maintain this software and are not responsible for its results.**

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Requirements

- Python 3.10+  
- Tkinter (usually included with Python)  

Optional (for packaging):  
- `pyinstaller`

---

## Installation

### 1. Create and activate a virtual environment (recommended)

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux / Mac
source .venv/bin/activate
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Running the app

```bash
python app.py
```

### 4. Packaging the app

```bash
pip install pyinstaller
pyinstaller --windowed --add-data "data;data" app.py
```


