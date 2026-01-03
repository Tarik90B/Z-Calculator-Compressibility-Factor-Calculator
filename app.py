"""
Z - Compressibility Factor Calculator
Copyright (c) 2026 Tarik

Independent implementation of AGA8 according to ISO 12213.
Licensed under the MIT License.
"""

import tkinter as tk
from frames.landing import LandingFrame
from frames.calculation import CalculationFrame
from frames.archive import ArchiveFrame
from tkinter import messagebox

def show_about():
    messagebox.showinfo(
        "About Z Calculator",
        "Z - Compressibility Factor Calculator\n"
        "Author: Tarik\n"
        "Version: 1.0\n\n"
        "License: MIT License\n"
        "Implements AGA8 (ISO 12213) equation\n"
        "ISO and AGA are not affiliated with this software"
    )

class App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Z - Compresibily factor calculator")
        self.geometry("1400x1000")

        menu = tk.Menu(self)
        self.config(menu=menu)

        menu.add_command(label="Home", command=lambda: self.show_frame("LandingFrame"))
        menu.add_command(label="Calculation", command=lambda: self.show_frame("CalculationFrame"))
        menu.add_command(label="Archive", command=lambda: self.show_frame("ArchiveFrame"))
        help_menu = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=show_about)


        container = tk.Frame(self)
        container.pack(fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
    
        self.frames = {}

        for FrameClass in (LandingFrame, CalculationFrame, ArchiveFrame):
            frame = FrameClass(container, self)
            self.frames[FrameClass.__name__] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("LandingFrame")

    def show_frame(self, frame_name):
        frame = self.frames[frame_name]
        frame.tkraise()

        if hasattr(frame, "refresh"):
            frame.refresh()

if __name__ == "__main__":
    app = App()
    app.mainloop()