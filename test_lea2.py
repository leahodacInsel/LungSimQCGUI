try:
    import tkinter as tk
    import tkinter.ttk as ttk
    from ttkthemes import ThemedTk

except ImportError:
    import Tkinter as tk
    import ttk

root = ThemedTk()
style = ttk.Style(root)

#    try also the 'clam' theme
style.theme_use('adapta')

common_bg = '#' + ''.join([hex(x)[2:].zfill(2) for x in (181, 26, 18)])  # RGB in dec
#    alternatively use the "more red" version of the common_bg as the indicatorcolor
#    sel_bg = '#' + ''.join([hex(x)[2:].zfill(2) for x in (221, 26, 18)])
common_fg = '#ffffff'  # pure white

rad_button = ttk.Radiobutton(root, text='abc')
rad_button.pack(expand=True, fill='both')

style_name = rad_button.winfo_class()
style.configure(style_name, foreground=common_fg, background=common_bg, indicatorcolor=common_bg)

style.map(style_name,
          foreground = [('disabled', common_fg),
                      ('pressed', common_fg),
                      ('active', common_fg)],
          background = [('disabled', common_bg),
                      ('pressed', '!focus', common_bg),
                      ('active', common_bg)],
          indicatorcolor=[('selected', common_bg),
                          ('pressed', common_bg)]

          )
root.mainloop()