from tkinter import Tk, filedialog

def get_save_filename() -> str:
    root = Tk()
    root.withdraw()

    filename = filedialog.asksaveasfilename(
        defaultextension=".png",
        filetypes=[("PNG image", "*.png"), ("TIF image", "*.tif")],
    )

    root.destroy()
    
    return filename