import tkinter as tk

class LandingFrame(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        button_frame = tk.Frame(self)
        button_frame.grid(row=0, column=0, sticky="nw", padx=0, pady=0)

        self.grid_columnconfigure(0, weight=1)
        tk.Label(self, text="Z-Factor Application", font=("Arial", 20, "bold"))\
            .grid(row=1, column=0, sticky="n", pady=20)