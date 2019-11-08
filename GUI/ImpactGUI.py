#!/usr/bin/env python3 

import sys; 
import pathlib;

src_path = pathlib.Path(__file__).parent.joinpath("src").absolute().as_posix()
if not src_path in sys.path:
    sys.path.append(src_path)
if not 'ImpactMainWindow' in sys.modules:
    ImpactMainWindow = __import__('ImpactMainWindow')
else:
    eval('import ImpactMainWindow')
    ImpactMainWindow = eval('reload(ImpactMainWindow)') 

root = ImpactMainWindow.ImpactMainWindow()
ImpactMainWindow.MyMenu(root)

root.update()
w  = root.winfo_width()
h  = root.winfo_height()
ws = root.winfo_screenwidth() # width of the screen
hs = root.winfo_screenheight() # height of the screen
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)

root.geometry('%dx%d+%d+%d' % (w, h, x, y))
root.resizable(width=True, height=True)
#root.protocol("WM_DELETE_WINDOW", quitConfirm)


root.mainloop()
root.console.stop()