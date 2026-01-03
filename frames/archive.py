import tkinter as tk
import csv
from tkinter import ttk
from utils import resource_path
from pathlib import Path

archive_file = resource_path("archive/z_factor_archive.csv", writable=True)
archive_file.parent.mkdir(exist_ok=True)
ARCHIVE_FILE = Path(archive_file)

class ArchiveFrame(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        button_frame = tk.Frame(self)
        button_frame.pack(fill="x", pady=5)

        self.controller = controller    
        self.tree = ttk.Treeview(self, show="headings")
        self.tree.pack(fill="both", expand=True)

        tk.Button(button_frame, text="Delete Selected Row", fg="red",
                  command=self.delete_selected_row).pack(side="left", padx=5)

    def refresh(self):
        self.load_archive()
      
    def load_archive(self):
        if not ARCHIVE_FILE.exists():
            return

        for item in self.tree.get_children():
            self.tree.delete(item)

        with ARCHIVE_FILE.open(newline="") as file:
            reader = csv.reader(file)
            headers = next(reader)

            self.tree["columns"] = headers
           
            first_col_width = 120   
            last_col_width = 100     
            middle_col_width = 50   

            for i, col in enumerate(headers):
                self.tree.heading(col, text=col)
        
                if i == 0: 
                    self.tree.column(col, width=first_col_width, stretch=False, anchor="w")
                elif i == len(headers) - 1: 
                    self.tree.column(col, width=last_col_width, stretch=False, anchor="center")
                else:  
                    self.tree.column(col, width=middle_col_width, stretch=False, anchor="center")

            for i, row in enumerate(reader):
                self.tree.insert("", "end", iid=str(i), values=row)

    def delete_selected_row(self):
        selected = self.tree.selection()
        if not selected:
            return  

        row_index = int(selected[0])
    
        self.tree.delete(selected[0])
  
        if ARCHIVE_FILE.exists():
            with ARCHIVE_FILE.open(newline="") as file:
                lines = list(csv.reader(file))

            if 0 <= row_index < len(lines) - 1:
                lines.pop(row_index + 1)

            with ARCHIVE_FILE.open("w", newline="") as file:
                writer = csv.writer(file)
                writer.writerows(lines)

        for new_iid, item in enumerate(self.tree.get_children()):
            self.tree.item(item, iid=str(new_iid))
   
        