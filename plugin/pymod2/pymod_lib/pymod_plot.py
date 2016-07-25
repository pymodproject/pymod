"""
Minimal Tkinter implementation of a graphical plotting system. Provides a limited sets of
functionalities with respect to Matplotlib, but since there are problems distributing and running
Matplotlib in many PyMOL builds, PyMod makes use of this plotting engine.
NumPy and Pmw are required.
"""

from Tkinter import *
import tkFileDialog
import tkMessageBox
import StringIO
import math
# Pmw.
global has_pmw
try:
    import Pmw
    has_pmw = True
except:
    has_pmw = False
# NumPy.
global has_numpy
try:
    import numpy
    has_numpy = True
except:
    has_numpy = False


###################################################################################################
# Plotting window.                                                                                #
###################################################################################################

# TODO: cambia il nome in Plot_window.
class Custom_plot_window(Toplevel):

    def __init__(self, parent, title = None, **configs):
        # Sets the window parameters.
        Toplevel.__init__(self, parent, **configs)
        self.geometry("850x600")
        self.minsize("540","381")

        if title:
            self.title(title)
        self.withdraw()

        # Top frame of the window containing: the title, the Canvas witht the plot and a column for the labels.
        self.plot_frame = Frame(self, bd=2,relief=GROOVE, bg="gray70")
        self.plot_frame.grid(row=0, column=0, sticky=E+W+S+N, padx=5, pady=5)
        self.plot_frame.columnconfigure(0, weight=1)
        self.plot_frame.rowconfigure(1, weight=1)

        # A frame which will contain the canvas with the plots.
        self.canvas_plot_frame = Frame(self.plot_frame)
        self.canvas_plot_frame.grid(row=1, column=0, sticky=E+W+S+N, padx=5, pady=5)

        self.plots_control_frames_list = []
        self.plots_control_frames_dict = {}

        # A bottom frame fo the window, containing some buttons to interact with the graph.
        self.control_frame = Frame(self)# , bg="gray")
        self.control_frame.grid(row=1, column=0, sticky=W)

        # The widgets will be packed in the 'show' method.
        self.controls_padding = 4
        self.controls_font = None
        self.controls_config = {"fg":"black", "padx" : self.controls_padding, "pady" : self.controls_padding, "font": self.controls_font}
        self.labels_pack_config = {"side":"left","pady":(0,5),"padx":(5,0)}
        self.buttons_pack_config = {"side":"left","pady":(0,5),"padx":(1,0)}

        # View options.
        self.view_label = Label(self.control_frame, text="View:", **self.controls_config)
        self.home_view_button = Button(self.control_frame, text="Home", command=self.return_to_home_view, **self.controls_config)
        self.data_view_button = Button(self.control_frame, text="Fit to Data", command=self.adapt_view_to_data_limits, **self.controls_config)
        # Zoom options.
        self.zoom_label = Label(self.control_frame, text="Zoom:", **self.controls_config)
        self.zoom_view_factor = 10
        self.zoom_out_button = Button(self.control_frame, text=" - ", command=self.zoom_out_view, **self.controls_config)
        self.zoom_in_button = Button(self.control_frame, text="+", command=self.zoom_in_view, **self.controls_config)
        # Move options.
        self.move_label = Label(self.control_frame, text="Move:", **self.controls_config)
        self.move_view_factor = 15
        self.move_left_button = Button(self.control_frame, text="Left", command=self.move_view_left, **self.controls_config)
        self.move_up_button = Button(self.control_frame, text="Up", command=self.move_view_up, **self.controls_config)
        self.move_down_button = Button(self.control_frame, text="Down", command=self.move_view_down, **self.controls_config)
        self.move_right_button = Button(self.control_frame, text="Right", command=self.move_view_right, **self.controls_config)

        # Interacting mode options.
        self.current_interaction_mode = StringVar()
        self.current_interaction_mode.set("interact") # initialize
        radiobutton_style_1 = {"foreground" : "black", "anchor" : "center", "font":self.controls_font,
                               "indicatoron": 0, "selectcolor" : "gray",
                               "padx" : self.controls_padding, "pady" : self.controls_padding}
        self.on_click_label = Label(self.control_frame, text="On click:", **self.controls_config)
        self.interact_mode_button = Radiobutton(self.control_frame, text="Interact", variable=self.current_interaction_mode, value="interact", command=self.activate_plot_interaction, **radiobutton_style_1)
        self.zoom_mode_button = Radiobutton(self.control_frame, text="Zoom", variable=self.current_interaction_mode, value="zoom", command=self.activate_plot_zooming, **radiobutton_style_1)

        # Export options.
        self.export_label = Label(self.control_frame, text="Export:", **self.controls_config)
        # The buttons will be packed in the 'show' method.
        self.export_to_csv_button = Button(self.control_frame, text="CSV", command=self.export_data, **self.controls_config)
        self.export_to_ps_button = Button(self.control_frame, text="PS", command=self.save_postscript, **self.controls_config)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)


    def build_plotting_area(self,
        use_plotting_controls = True,
        message_bar_initial_text = None,
        on_click_action=None,
        update_message_bar=False,
        message_bar_text_on_update="",
        message_bar_vars_on_update=(),
        x_label_text = "x",
        y_label_text = "y",
        hide_x_thicks = False,
        hide_y_thicks = False,
        hide_x_label = False,
        hide_y_label = False,
        highlight_points_on_click=True):
        """
        Builds the 'canvas_plot' which consists in the actual drawing area. It must be called after
        the initializing an object of this class.
        """
        # A label with the title of the plots canvas.
        if message_bar_initial_text:
            self.plot_title = Label(self.plot_frame, text=message_bar_initial_text, font="comic 10", bg="gray70")
            self.plot_title.grid(row=0, column=0, sticky=W, padx=5, pady=5, columnspan=2)

        self.canvas_plot = Canvas_plot(self.canvas_plot_frame,
                                       on_click_action=on_click_action,
                                       update_message_bar = update_message_bar,
                                       message_bar_text_on_update = message_bar_text_on_update,
                                       message_bar_vars_on_update = message_bar_vars_on_update,
                                       x_label_text = x_label_text,
                                       y_label_text = y_label_text,
                                       hide_x_thicks = hide_x_thicks,
                                       hide_y_thicks = hide_y_thicks,
                                       hide_x_label = hide_x_label,
                                       hide_y_label = hide_y_label,
                                       highlight_points_on_click = highlight_points_on_click,
                                       bg="gray50") # ,highlightthickness=0,insertbackground="red") # "gray40"
        self.canvas_plot.pack(expand=True, fill="both")

        # Builds a right column on the plot window. It will contain checkbuttons and legends to
        # show/hide each plot.
        self.use_plotting_controls = use_plotting_controls
        if self.use_plotting_controls:
            # A frame with the labels of the plots.
            self.labels_frame = None
            self.inner_labels_frame = None
            if has_pmw:
                # Use Pmw to build a frame with a scrollbar.
                self.labels_frame = Pmw.ScrolledFrame(self.plot_frame,
                                                      vscrollmode='dynamic',
                                                      hscrollmode='dynamic',
                                                      # horizflex="shrink",
                                                      usehullsize=1,
                                                      hull_width=200, # 250
                                                      hull_bg='white', clipper_bg='white'
                                                      )
                # This is the actual Frame where the content of the tab is going to be packed.
                self.inner_labels_frame = self.labels_frame.interior()
                self.labels_frame.grid(row=1, column=1, sticky=N+W+E+S, padx=(0,5), pady=5)
            else:
                # Builds the frame.
                self.labels_frame = Frame(self.plot_frame, bg="gray",bd=1,relief="groove")
                # Builds a scrollbar.
                self.labels_frame_scrollbar = Scrollbar(self.labels_frame)
                self.labels_frame_scrollbar.pack(side="right", fill="y")
                # Builds a canvas inside the labels frame.
                self.inner_labels_frame = Canvas(self.labels_frame, yscrollcommand=self.labels_frame_scrollbar.set)
                self.inner_labels_frame.pack(side="left",fill="both")
                # Configures the scrollbar.
                self.labels_frame_scrollbar.config(command=self.inner_labels_frame.yview)
                self.labels_frame.grid(row=1, column=1, sticky=N+W, padx=(0,5), pady=5)

            self.labels_config = {"bg":"white","bd":0,"font": "comic 9"}
            self.labels_pack_config = {"pady":0,"ipady":2}

            self.labels_title = Label(self.inner_labels_frame, text=" Plots List", anchor="w", **self.labels_config)
            self.labels_title.pack(fill="x",**self.labels_pack_config)


    def check_plotting_area(self):
        if not hasattr(self, "canvas_plot"):
            raise Exception("This 'Custom_plot_window' does not have a plotting area. Before starting to draw on it, use the 'build_plotting_area' method.")


    #################################################################
    # Add elements on the plotting area.                            #
    #################################################################

    def add_plot(self, x_data, y_data, **configs):
        """
        A wrapper method for the 'add_plot' method of the 'Canvas_plot' class.
        """
        self.check_plotting_area()
        if not configs.has_key("label"):
            configs["label"] = "a custom plot"
        plot = self.canvas_plot.add_plot(x_data, y_data, **configs)
        if self.use_plotting_controls:
            self.add_plot_control_frame(plot, configs["label"])


    def add_plot_control_frame(self, plot, label):
        """
        Builds a frame with a legend and controls for each plot.
        """
        self.check_plotting_area()
        # Build a plot control frame.
        plot_control_frame = Plot_control_frame(self.inner_labels_frame, bg="white")
        # Build its chekcbutton.
        cvar = IntVar()
        plot_checkbutton = Checkbutton(plot_control_frame, text=label, variable=cvar, highlightthickness=0,**self.labels_config)
        plot_checkbutton.var = cvar
        plot_checkbutton.configure(command = lambda x=plot_checkbutton: self.click_plot_checkbutton(x))
        # Add the checkbutton to the frame.
        plot_control_frame.add_checkbutton(plot_checkbutton)

        # Add a canvas with a line (or a dot) to work as a legend.
        legend_canvas_height, legend_canvas_width = plot_checkbutton.winfo_reqheight(),25
        legend_canvas = Canvas(plot_control_frame, width=legend_canvas_width, height=legend_canvas_height, bg="white", bd=0, highlightthickness=0,insertbackground="white")
        legend_canvas.create_line(0, float(legend_canvas_height)/2.0,legend_canvas_width,float(legend_canvas_height)/2.0, width=2, fill=plot.color)
        # Add a legend canvas.
        plot_control_frame.add_legend_canvas(legend_canvas)

        self.plots_control_frames_list.append(plot_control_frame)
        self.plots_control_frames_dict[plot_checkbutton] = plot


    def click_plot_checkbutton(self, checkbutton):
        """
        Shows or hides a plot in the plotting area when a user click on the plot checkbutton.
        """
        plot = self.plots_control_frames_dict[checkbutton]

        # Uncomment this to remove the last highlighted item when hiding or showing a plot.
        # self.canvas_plot.remove_last_clicked_item()

        if not checkbutton.var.get():
            # Hide the line of the plot and its points.
            for segment in plot.segments_list:
                self.canvas_plot.itemconfig(segment.id, state='hidden')
            map(lambda p: self.canvas_plot.itemconfig(p.id, state='hidden'), plot.get_points())
        else:
            # Show the line of the plot and its points.
            for segment in plot.segments_list:
                self.canvas_plot.itemconfig(segment.id, state='normal')
            map(lambda p: self.canvas_plot.itemconfig(p.id, state='normal'), plot.get_points())


    def draw_line(self, coords):
        """
        A wrapper method for the 'add_line' and 'show_line' methods of the 'Canvas_plot' class.
        """
        line = self.canvas_plot.add_line(coords)
        self.canvas_plot.show_line(line)


    def draw_label(self, label_text, label_coords):
        """
        A wrapper method for the 'add_label' and 'show_label' methods of the 'Canvas_plot' class.
        """
        label = self.canvas_plot.add_label(label_text, label_coords)
        self.canvas_plot.show_label(label)


    def show(self, **configs):
        """
        Called when the plot window containing all the plots must be showed.
        """
        self.check_plotting_area()
        self.canvas_plot.initialize(**configs)
        self.canvas_plot.show_all_objects()
        # Packs the plot control frames.
        for plot_control_frame in self.plots_control_frames_list:
            plot_control_frame.legend_canvas.pack(side="left", padx=5,**self.labels_pack_config)
            plot_control_frame.checkbutton.pack(side="left", **self.labels_pack_config)
            plot_control_frame.checkbutton.select()
            plot_control_frame.pack(expand=True, fill="both")
        # Packs the buttons.
        self.view_label.pack(**self.buttons_pack_config)
        self.home_view_button.pack(**self.buttons_pack_config)
        if self.canvas_plot.plots_list != []:
            self.data_view_button.pack(**self.buttons_pack_config)
        self.zoom_label.pack(**self.buttons_pack_config)
        self.zoom_out_button.pack(**self.buttons_pack_config)
        self.zoom_in_button.pack(**self.buttons_pack_config)
        self.move_label.pack(**self.buttons_pack_config)
        self.move_left_button.pack(**self.buttons_pack_config)
        self.move_up_button.pack(**self.buttons_pack_config)
        self.move_down_button.pack(**self.buttons_pack_config)
        self.move_right_button.pack(**self.buttons_pack_config)
        self.on_click_label.pack(**self.buttons_pack_config)
        self.interact_mode_button.pack(**self.buttons_pack_config)
        self.zoom_mode_button.pack(**self.buttons_pack_config)
        self.export_label.pack(**self.buttons_pack_config)
        if self.canvas_plot.plots_list != []:
            self.export_to_csv_button.pack(**self.buttons_pack_config)
        self.export_to_ps_button.pack(**self.buttons_pack_config)
        self.deiconify()


    #################################################################
    # Bottom frame controls events.                                 #
    #################################################################

    def return_to_home_view(self):
        self.canvas_plot.change_view(home_view=True)

    def adapt_view_to_data_limits(self):
        self.canvas_plot.change_view(x_limits="use_data_limits", y_limits="use_data_limits")

    def adapt_view_auto(self):
        self.canvas_plot.change_view(x_limits="auto", y_limits="auto")

    def activate_plot_interaction(self):
        self.canvas_plot.interaction_mode = self.current_interaction_mode.get()
        self.canvas_plot.configure(cursor="arrow")

    def activate_plot_zooming(self):
        self.canvas_plot.interaction_mode = self.current_interaction_mode.get()
        self.canvas_plot.configure(cursor="cross")

    def get_data_limits(self):
        return self.canvas_plot.current_plot_mins["x"], self.canvas_plot.current_plot_maxs["x"], self.canvas_plot.current_plot_mins["y"], self.canvas_plot.current_plot_maxs["y"]

    def zoom_in_view(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dx = (x_max - x_min)/self.zoom_view_factor
        dy = (y_max - y_min)/self.zoom_view_factor
        self.canvas_plot.change_view(x_limits=(x_min+dx, x_max-dx), y_limits=(y_min+dy, y_max-dy))

    def zoom_out_view(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dx = (x_max - x_min)/self.zoom_view_factor
        dy = (y_max - y_min)/self.zoom_view_factor
        self.canvas_plot.change_view(x_limits=(x_min-dx, x_max+dx), y_limits=(y_min-dy, y_max+dy))

    def move_view_left(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dx = (x_max - x_min)/self.move_view_factor
        self.canvas_plot.change_view(x_limits=(x_min-dx, x_max-dx), y_limits=(y_min, y_max))

    def move_view_up(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dy = (y_max - y_min)/self.move_view_factor
        self.canvas_plot.change_view(x_limits=(x_min, x_max), y_limits=(y_min+dy, y_max+dy))

    def move_view_down(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dy = (y_max - y_min)/self.move_view_factor
        self.canvas_plot.change_view(x_limits=(x_min, x_max), y_limits=(y_min-dy, y_max-dy))

    def move_view_right(self):
        x_min, x_max, y_min, y_max = self.get_data_limits()
        dx = (x_max - x_min)/self.move_view_factor
        self.canvas_plot.change_view(x_limits=(x_min+dx, x_max+dx), y_limits=(y_min, y_max))

    ##############################
    # Exporting images and data. #
    ##############################

    def export_data(self):
        if self.canvas_plot.plots_list != []:
            filepath = tkFileDialog.asksaveasfilename(filetypes=[("csv","*.csv")],parent=self)
            if not filepath == "":
                try:
                    output_file_handler = open(filepath,"w")
                    print >> output_file_handler, self.canvas_plot.export_to_csv()
                    output_file_handler.close()
                except:
                    pass
        else:
            tkMessageBox.showerror("Export Error", "No data to export.", parent=self)

    def save_postscript(self):
        filepath = tkFileDialog.asksaveasfilename(filetypes=[("ps","*.ps")],parent=self)
        if not filepath == "":
            try:
                # Lowers the stack level of the points in order to hide them in the .ps image.
                self.canvas_plot.tag_lower("type:point")
                self.canvas_plot.postscript(file=filepath, colormode='color')
                self.canvas_plot.set_canvas_stack()
            except:
                pass


class Plot_control_frame(Frame):
    """
    A class which will represent for each plot a frame containing a:
        - legend canvas, a canvas with a line with the same color as the respective line in the
          plot.
        - a checkbutton with the name of the line, if users click on it, the line (and its points)
          will be hidden in the plot.
    """
    def __init__(self, parent = None, **configs):
        Frame.__init__(self, parent, **configs)
        self.checkbutton = None

    def add_checkbutton(self, checkbutton):
        self.checkbutton = checkbutton

    def add_legend_canvas(self, legend_canvas):
        self.legend_canvas = legend_canvas


###################################################################################################
# Plotting canvas.                                                                                #
###################################################################################################

class Canvas_plot(Canvas):
    """
    A class to represent a plotting area in which plots, lines and points can be drawn.
    """

    plot_index = 0
    line_index = 0
    label_index = 0

    plot_colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "gray"]
    plot_color_index = 0

    def __init__(self, parent,
                 highlight_points_on_click, on_click_action,
                 update_message_bar, message_bar_text_on_update, message_bar_vars_on_update,
                 x_label_text, y_label_text,
                 hide_y_thicks, hide_x_thicks,
                 hide_x_label, hide_y_label,
                 **configs):
        Canvas.__init__(self, parent, **configs)

        self.plots_list = []
        self.plots_dict = {}

        self.lines_list = []
        self.lines_dict = {}

        self.labels_list = []
        self.labels_dict = {}

        # Interactions.
        self.on_click_action = on_click_action
        self.highlight_points_on_click = highlight_points_on_click
        self.update_message_bar = message_bar_text_on_update
        self.message_bar_text_on_update = message_bar_text_on_update
        self.message_bar_vars_on_update = message_bar_vars_on_update
        self.interaction_mode = "interact" # "zoom"

        # Plot appeareance.
        self.x_label_text = x_label_text
        self.y_label_text = y_label_text

        self.hide_x_thicks = hide_x_thicks
        self.hide_y_thicks = hide_y_thicks

        self.hide_x_label = hide_x_label
        self.hide_y_label = hide_y_label


    def change_plot_color_index(self):
        if self.plot_color_index == len(self.plot_colors) - 1:
            self.plot_color_index = 0
        else:
            self.plot_color_index += 1


    #################################################################
    # Initialization.                                               #
    #################################################################

    def initialize(self, x_limits="auto", y_limits="auto"):
        # Initializes a series of values needed to build the plot and all its components.
        self.axes_limits = {"x": x_limits, "y": y_limits} # "auto", "use_data_limits", (min, max)
        self.set_axes_input_limits()
        self.draw_extreme_thicks_lines = False

        self.min_data_values = {"x":None,"y":None}
        self.max_data_values = {"x":None,"y":None}
        # self.linspaces = {"x":None,"y":None}
        self.thicks_spaces = {"x":None,"y":None}
        self.thicks_labels_spaces = {"x":None,"y":None}
        self.original_plot_mins = {"x":None,"y":None}
        self.original_plot_maxs = {"x":None,"y":None}
        self.current_plot_mins = {"x":None,"y":None}
        self.current_plot_maxs = {"x":None,"y":None}
        self.labels = {"x":None,"y":None}
        self.numbers_of_thicks = {"x":None,"y":None}

        self.axes_labels_x_coords = {"x":None,"y":None}
        self.axes_labels_y_coords = {"x":None,"y":None}
        self.thicks_labels_x_coords = {"x":None,"y":None}
        self.thicks_labels_y_coords = {"x":None,"y":None}
        self.thicks_lines_x_coords = {"x":None,"y":None}
        self.thicks_lines_y_coords = {"x":None,"y":None}

        self.axes_labels_texts = {"x":self.x_label_text,"y":self.y_label_text}
        self.hide_axes_thicks = {"x": self.hide_x_thicks, "y":self.hide_y_thicks}
        self.hide_axes_labels = {"x": self.hide_x_label, "y":self.hide_y_label}

        # Initiliazes a series of coordinates which will be used to display the elements of the plot.
        self.set_canvas_plot_coords(initialize=True)

        ###########################################
        # Builds the white-colored plotting area. #
        ###########################################
        self.plotting_area = self.create_rectangle(self.plotting_area_x0,self.plotting_area_y0,self.plotting_area_x1,self.plotting_area_y1, fill="white", tags=("plotting_area"))
        self.tag_bind(self.plotting_area, '<Button-1>', self.click_on_plotting_area)
        # Build a "frame" around the plotting area, so that the objects plotted outside it will be
        # hidden.
        self.plotting_area_frame_pieces = [None,None,None,None]
        bg_color = self["background"] # self["background"]
        self.plotting_area_frame_pieces[0] = self.create_rectangle(0,0,self.width,self.plotting_area_dy, fill=bg_color, outline=bg_color, tags=("plotting_area_frame_piece"))
        self.plotting_area_frame_pieces[1] = self.create_rectangle(self.width-self.plotting_area_dx,0,self.width,self.height, fill=bg_color, outline=bg_color, tags=("plotting_area_frame_piece"))
        self.plotting_area_frame_pieces[2] = self.create_rectangle(0,self.height-self.plotting_area_dy,self.width,self.height, fill=bg_color, outline=bg_color, tags=("plotting_area_frame_piece"))
        self.plotting_area_frame_pieces[3] = self.create_rectangle(0,0,self.plotting_area_dx,self.height, fill=bg_color, outline=bg_color, tags=("plotting_area_frame_piece"))
        # Borders of the plotting area.
        self.plotting_area_borders = [None,None,None,None]
        bd_color = "black"
        self.plotting_area_frame_pieces[0] = self.create_line(self.plotting_area_dx,self.plotting_area_dy,self.width-self.plotting_area_dx,self.plotting_area_dy, fill=bd_color, tags=("plotting_area_border"))
        self.plotting_area_frame_pieces[1] = self.create_line(self.width-self.plotting_area_dx,self.plotting_area_dy,self.width-self.plotting_area_dx,self.height-self.plotting_area_dy, fill=bd_color, tags=("plotting_area_border"))
        self.plotting_area_frame_pieces[2] = self.create_line(self.plotting_area_dx,self.height-self.plotting_area_dy,self.width-self.plotting_area_dx,self.height-self.plotting_area_dy, fill=bd_color, tags=("plotting_area_border"))
        self.plotting_area_frame_pieces[3] = self.create_line(self.plotting_area_dx,self.plotting_area_dy,self.plotting_area_dx,self.height-self.plotting_area_dy, fill=bd_color, tags=("plotting_area_border"))

        #####################################
        # Allows the zooming functionality. #
        #####################################
        self.zx = self.zy = 0
        self.tag_bind(self.plotting_area, "<ButtonRelease-1>", self.on_mouse_release)
        self.tag_bind(self.plotting_area, "<B1-Motion>", self.on_mouse_drag)
        self.zlines = [None,None,None,None]
        self.zlines[0] = self.create_line(0,0,0,0) # Top horizontal line.
        self.zlines[1] = self.create_line(0,0,0,0) # Right vertical line.
        self.zlines[2] = self.create_line(0,0,0,0) # Bottom horizontal line.
        self.zlines[3] = self.create_line(0,0,0,0) # Left vertical line.
        for zline in self.zlines:
            self.itemconfig(zline, state='hidden')
        self.showing_zlines = False

        ###################
        # Build the axes. #
        ###################
        self.build_axis("x", mode="initialize")
        self.build_axis("y", mode="initialize")

        # Configure the plot to be resizable.
        self.bind("<Configure>", self.on_resize)

        # Configure plot interactions.
        self.last_clicked_item = None


    def set_axes_input_limits(self):
        self.input_min_limits = {"x":None, "y":None}
        self.input_max_limits = {"x":None, "y":None}
        if hasattr(self.axes_limits["x"], "__iter__"):
            self.input_min_limits["x"] = self.axes_limits["x"][0]
            self.input_max_limits["x"] = self.axes_limits["x"][1]
        if hasattr(self.axes_limits["y"], "__iter__"):
            self.input_min_limits["y"] = self.axes_limits["y"][0]
            self.input_max_limits["y"] = self.axes_limits["y"][1]


    def change_view(self, x_limits="auto", y_limits="auto", home_view=False):
        if home_view:
            x_limits = self.original_axes_limits["x"]
            y_limits = self.original_axes_limits["y"]
        # Change the axis range and thicks, that is, rebuilds the axes.
        self.axes_limits = {"x": x_limits, "y": y_limits}
        self.set_axes_input_limits()
        self.build_axis("x", mode="rebuild")
        self.build_axis("y", mode="rebuild")
        # Replots the objects on the plot in order to adapt them to the new view.
        self.remove_last_clicked_item()
        self.show_all_objects()


    def set_canvas_plot_coords(self, initialize=False):
        # Gets the current dimensions of the Canvas widget.
        self.width = self.winfo_reqwidth()
        self.height = self.winfo_reqheight()
        # Relative margin between the Canvas border and the inner white-colored drawing area.
        self.relative_x_margin = 0.12
        self.relative_y_margin = 0.10
        self.plotting_area_dx = self.width*self.relative_x_margin
        self.plotting_area_dy = self.height*self.relative_y_margin
        # Coordinates to build the plotting area.
        self.plotting_area_x0 = self.plotting_area_dx
        self.plotting_area_y0 = self.plotting_area_dy
        self.plotting_area_x1 = self.width-self.plotting_area_dx
        self.plotting_area_y1 = self.height-self.plotting_area_dy

        # X labels coords.
        self.axes_labels_x_coords["x"] = self.width/2 # [(self.width)/2, (self.height-self.plotting_area_dy)+self.x_label_offset]
        self.axes_labels_y_coords["x"] = self.plotting_area_y1+25
        # Y labels coords.
        if initialize:
            self.y_label_offset = 55
            self.y_thicks_offset = 20
        self.axes_labels_x_coords["y"] = self.plotting_area_dx-self.y_label_offset
        self.axes_labels_y_coords["y"] = self.height/2
        # X thicks.
        self.thicks_labels_y_coords["x"] = self.plotting_area_y1+10
        self.thicks_lines_y_coords["x"] = self.plotting_area_y1
        # Y thicks.
        self.thicks_labels_x_coords["y"] = self.plotting_area_dx-self.y_thicks_offset
        self.thicks_lines_x_coords["y"] = self.plotting_area_dx


    #################################################################
    # Builds and displays axes.                                     #
    #################################################################

    def set_limit_data_values(self, axis):
        """
        Gets the min and max values from the data stored in the 'Custom_plot' objects added through
        the 'add_plot' method.
        """
        if self.plots_list != []:
            all_axis_data = reduce(lambda v1,v2: v1+v2, [plot.get_data_coords(component=axis) for plot in self.plots_list])
            if axis == "y":
                self.min_data_values[axis] = min(filter(lambda v: v != None, all_axis_data))
                self.max_data_values[axis] = max(filter(lambda v: v != None, all_axis_data))
            else:
                self.min_data_values[axis] = min(all_axis_data)
                self.max_data_values[axis] = max(all_axis_data)


    def build_axis(self, axis, mode="initialize"):
        """
        Builds a series of labels and thicks for an axis.
        When the axes are built for the first time, the 'change_view' attribute has to be set to
        'False'.
        """
        if mode in ("initialize", "update"):
            self.set_limit_data_values(axis)

        # Builds lists of evenly spaced thicks and labels and gets the current limits for the
        # plotting area.
        if (self.axes_limits[axis] == "auto" or self.axes_limits[axis] == "use_data_limits") and self.plots_list == []:
            raise Exception("The plotting area limits can not be defined using min and max values from plots data since no plot has been added through the 'add_plot' method.")
        if self.axes_limits[axis] == "auto":
            self.thicks_spaces[axis], self.thicks_labels_spaces[axis] = get_plotting_interval(self.min_data_values[axis], self.max_data_values[axis], mode="expand")
            self.current_plot_mins[axis], self.current_plot_maxs[axis] = self.thicks_spaces[axis][0], self.thicks_spaces[axis][-1]
        elif self.axes_limits[axis] == "use_data_limits":
            self.thicks_spaces[axis], self.thicks_labels_spaces[axis] = get_plotting_interval(self.min_data_values[axis], self.max_data_values[axis], mode="use_min_and_max")
            self.current_plot_mins[axis], self.current_plot_maxs[axis] = self.thicks_spaces[axis][0], self.thicks_spaces[axis][-1]
        # Use limits specified in the 'initialize' method of this class.
        elif hasattr(self.axes_limits[axis], "__iter__"):
            self.thicks_spaces[axis], self.thicks_labels_spaces[axis] = get_plotting_interval(self.input_min_limits[axis], self.input_max_limits[axis], mode="use_min_and_max")
            self.current_plot_mins[axis], self.current_plot_maxs[axis] = self.thicks_spaces[axis][0], self.thicks_spaces[axis][-1]

        if mode == "initialize":
            self.original_plot_mins[axis] = self.current_plot_mins[axis]
            self.original_plot_maxs[axis] = self.current_plot_maxs[axis]
            self.original_axes_limits = self.axes_limits.copy()

        if self.current_plot_maxs[axis] == self.current_plot_mins[axis]:
            self.current_plot_maxs[axis] *= 1.05
            if self.current_plot_maxs[axis] == 0:
                self.current_plot_maxs[axis] += 0.05

        # Builds the axis label.
        if mode in ("initialize", "update"):
            if axis == "y":
                text_to_display = "\n".join(self.axes_labels_texts[axis])
            else:
                text_to_display = self.axes_labels_texts[axis]
            self.labels[axis] = self.create_text(self.axes_labels_x_coords[axis],self.axes_labels_y_coords[axis], text=text_to_display,tags = "%s_label" % (axis))
            if self.hide_axes_labels[axis]:
                self.itemconfig(self.labels[axis], state='hidden')

        # Thicks and their labels.
        # Remove old thicks and labels if changing view or updating.
        if mode in ("rebuild", "update"):
            for thick_label in self.find_withtag("%s_thick_label" % (axis)):
                self.delete(thick_label)
            for thick_line in self.find_withtag("%s_thick_line" % (axis)):
                self.delete(thick_line)

        # Builds new thicks and labels.
        self.numbers_of_thicks[axis] = len(self.thicks_spaces[axis]) - 1
        thick_counter = 0
        for thick_text, thick_variable_coord in zip(self.thicks_labels_spaces[axis], self.get_thick_variable_coords(axis)):
            # Builds the label of the thick.
            thick_label_x, thick_label_y = self.build_thick_label_coords(axis, thick_variable_coord)
            thick_label = self.create_text(thick_label_x, thick_label_y, text=thick_text, tags="%s_thick_label"%(axis))
            if self.hide_axes_thicks[axis]:
                self.itemconfig(thick_label, state='hidden')
            # Builds the line of the thick.
            thick_line_x0, thick_line_y0, thick_line_x1, thick_line_y1 = self.build_thick_line_coords(axis, thick_variable_coord)
            thick_line = self.create_line(thick_line_x0, thick_line_y0, thick_line_x1, thick_line_y1, tags="%s_thick_line"%(axis))
            if self.hide_axes_thicks[axis]:
                self.itemconfig(thick_line, state='hidden')
            if not self.draw_extreme_thicks_lines:
                if thick_counter == 0 or thick_counter == self.numbers_of_thicks[axis]:
                    self.itemconfig(thick_line, state='hidden')
            thick_counter += 1


    def get_thick_variable_coords(self, axis):
        return [self.convert_single_coordinate(axis, t) for t in self.thicks_spaces[axis]]

    def build_thick_label_coords(self, axis, variable_coord):
        if axis == "x":
            return variable_coord, self.thicks_labels_y_coords["x"]
        elif axis == "y":
            return self.thicks_labels_x_coords["y"], variable_coord

    def build_thick_line_coords(self, axis, variable_coord):
        if axis == "x":
            return variable_coord, self.thicks_lines_y_coords["x"], variable_coord, self.thicks_lines_y_coords["x"]-5
        elif axis == "y":
            return self.thicks_lines_x_coords["y"], variable_coord, self.thicks_lines_x_coords["y"]+5, variable_coord


    #################################################################
    # Resizing events.                                              #
    #################################################################

    def on_resize(self, event):
        # Determine the ratio of old width/height to new width/height.
        wscale = float(event.width)/self.width
        hscale = float(event.height)/self.height
        wdiff = float(event.width) - self.width
        hdiff = float(event.height) - self.height
        # Resize the plotting canvas.
        self.config(width=event.width, height=event.height)
        # Rescale the plotting area and all the objects tagged with the "resize" tag.
        self.scale("plotting_area",0,0,wscale,hscale)
        self.scale("plotting_area_frame_piece",0,0,wscale,hscale)
        self.scale("plotting_area_border",0,0,wscale,hscale)
        self.scale("resize",0,0,wscale,hscale)

        # Resacale other objects.
        self.set_canvas_plot_coords()
        # Repositions the x and y axes labels.
        self.coords(self.labels["x"],self.axes_labels_x_coords["x"], self.axes_labels_y_coords["x"])
        self.coords(self.labels["y"],self.axes_labels_x_coords["y"], self.axes_labels_y_coords["y"])

        # Repositions thicks.
        for axis in ("x","y"):
            thicks_labels = self.find_withtag("%s_thick_label" % (axis))
            thicks_lines = self.find_withtag("%s_thick_line" % (axis))
            for thick_label_id, thick_line_id, thick_variable_coord in zip(thicks_labels, thicks_lines, self.get_thick_variable_coords(axis)):
                # Repositions the labels of the thicks.
                thick_label_x, thick_label_y = self.build_thick_label_coords(axis,thick_variable_coord)
                self.coords(thick_label_id, thick_label_x, thick_label_y)
                # Repositions the lines of the thicks.
                thick_line_x0, thick_line_y0, thick_line_x1, thick_line_y1 = self.build_thick_line_coords(axis, thick_variable_coord)
                self.coords(thick_line_id, thick_line_x0, thick_line_y0, thick_line_x1, thick_line_y1)

        # Resize highlighted point.
        highlighted = self.find_withtag("highlighted")
        for h in highlighted:
            cx,cy = self.get_point_center(self.coords(h))
            dx = ((cx-0)*wscale-0) - cx
            dy = ((cy-0)*hscale+0) - cy
            self.move(h, dx, dy)

        # Finally, updates the dimension according to the event.
        self.width = event.width
        self.height = event.height


    #################################################################
    # Add new object to the plotting area.                          #
    #################################################################

    def add_plot(self, x_data, y_data, **configs):
        """
        Adds a plot to be drawn. The canvas coordinates of its components wil be initially set to 0.
        They will be changed when the 'show_plot' method is called.
        """
        # Prepares x data.
        if x_data == None:
            x_data = range(0,len(y_data))
        # Prepares the additional data.
        if not configs.has_key("additional_data") or (configs.has_key("additional_data") and configs["additional_data"] == None):
            configs["additional_data"] = [{None:None}]*len(y_data)

        # Adds the new plot.
        if not configs.has_key("color"):
            configs["color"] = self.plot_colors[self.plot_color_index]
            self.change_plot_color_index()
        new_plot = Custom_plot(x_data, y_data, plot_id=self.plot_index, **configs)
        self.plots_list.append(new_plot)
        self.plots_dict[self.plot_index] = new_plot
        self.plot_index += 1

        # Builds the line with all its segments (that is, series of points separated by points with
        # 'None' y values).
        for segment in new_plot.segments_list:
            # TODO: remove.
            # line_coords = []
            # for point in segment.points:
            #     line_coords.extend([0,0])
            # # Segments with only one points.
            # if len(line_coords) < 4:
            #     line_coords.extend([0, 0])
            pid = self.create_line([0,0,0,0])
            segment.id = pid
            self.itemconfig(pid, width=2, fill=new_plot.color, tags=("plot","resize"))
            self.tag_bind(pid, '<Button-1>', self.click_on_plot)
            self.tag_bind(pid, "<ButtonRelease-1>", self.on_mouse_release)
            self.tag_bind(pid, "<B1-Motion>", self.on_mouse_drag)

        # Builds the points.
        for point in new_plot.get_points(exclude_none_values=True):
            point.id = self.create_circle(0, 0, radius=2, state="normal",
                                          outline="black",
                                          width=0, # activewidth=5, activefill="red", activeoutline="blue", fill=plot.color, outline="black")
                                          )
            new_plot.points_dict[point.id] = point
            self.itemconfig(point.id, tags=("type:point", "plot:%s" % (new_plot.id), "point:%s" % (point.id), "resize"))
            self.tag_bind(point.id, "<Button-1>", self.click_on_point)
            self.tag_bind(point.id, "<ButtonRelease-1>", self.on_mouse_release)
            self.tag_bind(point.id, "<B1-Motion>", self.on_mouse_drag)

        return new_plot


    #################################################################
    # Drawing.                                                      #
    #################################################################

    def show_all_objects(self):
        for plot in self.plots_list:
            self.show_plot(plot)
        for line in self.lines_list:
            self.show_line(line)
        for label in self.labels_list:
            self.show_label(label)
        self.set_canvas_stack()


    def set_canvas_stack(self):
        """
        # Set the "z-level" of the items in the plot in order for it to be displayed correctly.
        """
        self.tag_raise("plot")
        self.tag_raise("type:point")
        self.tag_raise("plotting_area_frame_piece")
        self.tag_raise("plotting_area_frame_piece")
        self.tag_raise(self.labels["x"])
        self.tag_raise(self.labels["y"])
        self.tag_raise(self.axes_labels_x_coords["x"])
        self.tag_raise(self.axes_labels_x_coords["y"])
        self.tag_raise("x_thick_label")
        self.tag_raise("y_thick_label")
        self.tag_raise("label")
        self.tag_raise("plotting_area_border")


    def show_plot(self, plot):
        """
        Actually draws a plot.
        """
        # Assign to each data point of the plot the corresponding coordinates in the 'canvas'
        # system.
        for point in plot.get_points():
            point.xc, point.yc = self.convert_coordinates(point.xd, point.yd)

        # Builds the line with all its segments (that is, series of points separated by points with
        # 'None' y values).
        for segment in plot.segments_list:
            line_coords = []
            for point in segment.points:
                line_coords.extend([point.xc,point.yc])
            # Segments with only one points.
            if len(line_coords) < 4:
                line_coords.extend([segment.points[0].xc+1, segment.points[0].yc+1])
            self.coords(segment.id, *line_coords)
            self.tag_raise(segment.id)

        # Builds the points.
        for point in plot.get_points(exclude_none_values=True):
            coords = self.get_circle_coords(point.xc, point.yc, radius=2)
            self.coords(point.id, *coords)
            self.tag_raise(point.id)


    def add_line(self, line_coords):
        # Stores information about the line.
        line = Custom_line(line_coords)
        line.id = self.create_line([0,0,0,0], tags=("line","resize"))
        self.lines_list.append(line)
        self.lines_dict[self.line_index] = line
        self.line_index += 1
        return line

    def show_line(self, line):
        """
        Draws a line on the plot. The 'line_coords' are supplied in the data coordinates systems and
        will be converted in the canvas coordinate system.
        """
        # Gets the canvas coordinates of the line.
        line_coords_c = [None,None,None,None]
        line_coords_c[0], line_coords_c[1] = self.convert_coordinates(line.coordinates[0], line.coordinates[1])
        line_coords_c[2], line_coords_c[3] = self.convert_coordinates(line.coordinates[2], line.coordinates[3])
        self.coords(line.id, *line_coords_c)
        self.tag_raise(line.id)


    def add_label(self, label_text, label_coords):
        # Builds the label.
        label = Custom_label(label_text, label_coords)
        label.id = self.create_text(0, 0, anchor="w", text=label_text,tags = ("label", "resize"))
        # Stores information about the label.
        self.labels_list.append(label)
        self.labels_dict[self.label_index] = label
        self.label_index += 1
        return label

    def show_label(self, label):
        """
        Draws a label on the plot.
        """
        # Gets the canvas coordinates of the label.
        label_coords_c = [None,None]
        label_coords_c[0], label_coords_c[1] = self.convert_coordinates(label.coordinates[0], label.coordinates[1])
        self.coords(label.id, *label_coords_c)
        self.tag_raise(label.id)


    #############################################
    # Converts coordinates from  the 'data'     #
    # coordinate system to the 'canvas' one.    #
    #############################################

    def convert_coordinates(self, xd, yd):
        """
        Converts coordinates from  the 'data' coordinate system to the 'canvas' coordinate system,
        so that the data can be displayed correctly in the drawing area.
        """
        if xd != None and yd != None:
            coords = self.coords(self.plotting_area) # Returns something like: [38.0, 26.700000000000003, 342.0, 240.3]
            return self.convert_x_coordinate(xd, coords), self.convert_y_coordinate(yd, coords)
        else:
            return None, None

    def convert_single_coordinate(self, axis, data1, data2=0):
        """
        Similar to the 'convert_coordinates' method, but this also allows to convert coordinates
        from one axis only.
        """
        coords = self.coords(self.plotting_area)
        if axis == "x":
            if data1 != None:
                return self.convert_x_coordinate(data1, coords)
        elif axis == "y":
            if data1 != None:
                return self.convert_y_coordinate(data1, coords)
        elif axis == "xy":
            if data1 != None and data2 != None:
                return self.convert_x_coordinate(data1, coords), self.convert_y_coordinate(data2, coords)

    def convert_x_coordinate(self, xd, coords):
        xd = float(xd)
        dxc = coords[0]
        lc = coords[2] - coords[0]
        ld = self.current_plot_maxs["x"] - self.current_plot_mins["x"]
        xc = dxc + ((xd-self.current_plot_mins["x"])/ld) * lc # x 'canvas' coordinate.
        return xc

    def convert_y_coordinate(self, yd, coords):
        yd = float(yd)
        dyc = coords[1] # Top left y 'canvas' coordinate of the drawing area.
        hc = coords[3] - coords[1] # Gets the height of the plotting area, in 'canvas' units.
        hd = self.current_plot_maxs["y"] - self.current_plot_mins["y"]
        yc = dyc + hc - ((yd-self.current_plot_mins["y"])/hd) * hc # y 'canvas' coordinate.
        return yc

    #############################################
    # Converts coordinates from  the 'canvas'   #
    # coordinate system to the 'data' one.      #
    #############################################

    def convert_x_cd(self, xc, coords): # 'cd' stands for 'canvas to data'.
        xc = float(xc)
        dxc = coords[0]
        lc = coords[2] - coords[0]
        ld = self.current_plot_maxs["x"] - self.current_plot_mins["x"]
        xd = (xc - dxc)*(ld/lc) + self.current_plot_mins["x"] # x 'data' coordinate.
        return xd

    def convert_y_cd(self, yc, coords):
        yc = float(yc)
        dyc = coords[1] # Top left y 'canvas' coordinate of the drawing area.
        hc = coords[3] - coords[1] # Gets the height of the plotting area, in 'canvas' units.
        hd = self.current_plot_maxs["y"] - self.current_plot_mins["y"]
        yd = (dyc + hc - yc)*(hd/hc) + self.current_plot_mins["y"] # y 'data' coordinate.
        return yd

    #################
    # Draw circles. #
    #################

    def create_circle(self, xc, yc, radius=3, **configs):
        radius = float(radius)
        coords = self.get_circle_coords(xc, yc, radius)
        return self.create_oval(coords[0],coords[1],coords[2],coords[3], **configs)

    def get_circle_coords(self, xc, yc, radius=3):
        return xc-radius,yc-radius,xc+radius,yc+radius


    #################################################################
    # Events and interactions.                                      #
    #################################################################

    def click_on_plotting_area(self, event):
        if self.interaction_mode == "interact":
            items = self.find_items(event)
            nearest_point_item = self.get_nearest_point(event, items)
            self.remove_last_clicked_item()
            if nearest_point_item:
                self.activate_point(event, nearest_point_item)
        elif self.interaction_mode == "zoom":
            self.zx = event.x
            self.zy = event.y

    def click_on_plot(self, event):
        if self.interaction_mode == "interact":
            plot = self.find_closest(self.canvasx(event.x), self.canvasy(event.y))[0]
            self.remove_last_clicked_item()
            if not self.plots_dict.has_key(plot):
                self.click_on_plotting_area(event)
            else:
                plot_object = self.plots_dict[plot]
                nearest_point_item = self.get_nearest_point(event, plot_object.points_dict.keys(), distance_along="x")
                if nearest_point_item:
                    self.activate_point(event, nearest_point_item)
        elif self.interaction_mode == "zoom":
            self.zx = event.x
            self.zy = event.y

    def click_on_point(self, event):
        if self.interaction_mode == "interact":
            items = self.find_items(event)
            nearest_point_item = self.get_nearest_point(event, items)
            self.remove_last_clicked_item()
            if nearest_point_item:
                self.activate_point(event, nearest_point_item)
        elif self.interaction_mode == "zoom":
            self.zx = event.x
            self.zy = event.y


    #########################################
    # Methods to build a zooming rectangle. #
    #########################################
    def adjust_coordinates_to_area(self, x0, y0, x1, y1):
        if self.plotting_area_x1 <= x1:
            x1 = self.plotting_area_x1
        if self.plotting_area_x0 >= x1:
            x1 = self.plotting_area_x0
        if self.plotting_area_y1 <= y1:
            y1 = self.plotting_area_y1
        if self.plotting_area_y0 >= y1:
            y1 = self.plotting_area_y0
        return x0, y0, x1, y1

    def order_coordinates(self, c0, c1):
        if c0 > c1:
            return c1,c0
        else:
            return c0,c1

    def on_mouse_drag(self, event):
        if self.interaction_mode == "zoom":
            if not self.showing_zlines:
                for zline in self.zlines:
                    self.itemconfig(zline, state='normal')
                    self.tag_raise(zline)
                self.showing_zlines = True
            x0,y0 = (self.zx, self.zy)
            x1,y1 = (event.x, event.y)
            x0, y0, x1, y1 = self.adjust_coordinates_to_area(x0, y0, x1, y1)
            self.coords(self.zlines[0], x0,y0,x1,y0)
            self.coords(self.zlines[1], x1,y0,x1,y1)
            self.coords(self.zlines[2], x0,y1,x1,y1)
            self.coords(self.zlines[3], x0,y0,x0,y1)

    def on_mouse_release(self, event):
        if self.interaction_mode == "zoom":
            x0,y0 = (self.zx, self.zy)
            x1,y1 = (event.x, event.y)
            # x0, y0, x1, y1 = self.adjust_coordinates_to_area(x0, y0, x1, y1)
            for zline in self.zlines:
                self.itemconfig(zline, state='hidden')
            self.showing_zlines = False
            # Converts the canvas coordinates to data coordinates.
            coords = self.coords(self.plotting_area)
            xd0, xd1 = self.order_coordinates(self.convert_x_cd(x0, coords),self.convert_x_cd(x1, coords))
            yd0, yd1 = self.order_coordinates(self.convert_y_cd(y0, coords), self.convert_y_cd(y1, coords))
            if xd0 <= self.current_plot_mins["x"]:
                xd0 = self.current_plot_mins["x"]
            if xd1 >= self.current_plot_maxs["x"]:
                xd1 = self.current_plot_maxs["x"]
            if yd0 <= self.current_plot_mins["y"]:
                yd0 = self.current_plot_mins["y"]
            if yd1 >= self.current_plot_maxs["y"]:
                yd1 = self.current_plot_maxs["y"]
            # Checks if the rectangle drawn the user is large enough, before zooming.
            if xd1 - xd0 < (self.current_plot_maxs["x"] - self.current_plot_mins["x"])/100:
                return False
            if yd1 - yd0 < (self.current_plot_maxs["y"] - self.current_plot_mins["y"])/100:
                return False
            # Actually zooms the plot area using the new min and max defined using the zooming
            # rectangle.
            self.change_view(x_limits=(xd0,xd1), y_limits=(yd0,yd1))


    ############################################
    # Methods to find items on the canvas plot #
    # near the clicked region.                 #
    ############################################

    def get_event_coords(self, event):
        return self.canvasx(event.x), self.canvasy(event.y)


    def find_items(self, event, search_mode="overlapping", find_range=5):
        xcf, ycf = self.get_event_coords(event)
        items_found = None
        if search_mode == "overlapping":
            items_found = self.find_overlapping(xcf-find_range, ycf-find_range,xcf+find_range, ycf+find_range)
        elif search_mode == "enclosed":
            items_found = self.find_enclosed(xcf-find_range, ycf-find_range,xcf+find_range, ycf+find_range)
        return items_found


    def get_nearest_point(self, event, items, distance_along="xy"):
        points_list = filter(lambda i: self.is_point(i), items)
        distances_list = []
        xcf, ycf = self.get_event_coords(event)
        for point in points_list:
            # Gets the 'canvas' coordinates of center of the point.
            coords = self.coords(point)
            pxc = coords[0] + (coords[2] - coords[0])/2.0
            pyc = coords[1] + (coords[3] - coords[1])/2.0
            # Finds the distance from the center of the point to the point on the canvas where the
            # click occurred.
            distance = None
            if distance_along == "xy":
                distance = math.sqrt(pow(xcf-pxc,2)+pow(ycf-pyc,2))
            elif distance_along == "x":
                distance = abs(xcf-pxc)
            elif distance_along == "y":
                distance = abs(ycf-pyc)
            distances_list.append(distance)

        if points_list:
            min_distance = min(distances_list)
            for point, distance in zip(points_list, distances_list):
                p = self.get_point_from_item(point)
                if distance == min_distance:
                    return point
        else:
            return None


    def is_point(self, item):
        for tag in self.gettags(item):
            if tag[5:] =="point":
                return True


    def get_plot_window_message_bar(self):
        canvas_plot_frame = self.nametowidget(self.winfo_parent())
        plot_frame = self.nametowidget(canvas_plot_frame.winfo_parent())
        root = self.nametowidget(plot_frame.winfo_parent())
        return root.plot_title


    def activate_point(self, event, point_item):
        point, plot = self.get_point_from_item(point_item, get_plot=True)

        # Launch an event.
        if hasattr(self.on_click_action, "__call__"):
            self.on_click_action(point, plot)

        # Updates the message bar.
        if self.update_message_bar:
            message_bar = self.get_plot_window_message_bar()
            message_bar_update_tuple = tuple([point.additional_data[additional_data_key] for additional_data_key in self.message_bar_vars_on_update])
            updated_text = self.message_bar_text_on_update.replace("__plot_name__", plot.label)
            updated_text = updated_text.replace("__x__",str(point.xd))
            updated_text = updated_text.replace("__y__",str(point.yd))
            message_bar["text"] = updated_text % message_bar_update_tuple

        # Highlights point being clicked in the ploting area.
        if self.highlight_points_on_click:
            coords = self.coords(point.id)
            x0, y0, x1, y1 = self.get_highlighted_point_coords(coords)
            self.last_clicked_item = self.create_oval(x0, y0, x1, y1, state="normal",
                                                      outline="black", width=1, fill="yellow",tag="highlighted")

    def get_highlighted_point_coords(self, point_coords):
        r = 3
        cx, cy = self.get_point_center(point_coords)
        x0 = cx-r
        x1 = cx+r
        y0 = cy-r
        y1 = cy+r
        return x0, y0, x1, y1

    def get_point_center(self, point_coords):
        cx = (point_coords[2] + point_coords[0])/2.0
        cy = (point_coords[3] + point_coords[1])/2.0
        return cx, cy

    def remove_last_clicked_item(self):
        if self.last_clicked_item:
            self.delete(self.last_clicked_item)


    def get_point_from_item(self, item, get_plot=False):
        line = None
        point = None
        for tag in self.gettags(item):
            if tag[0:4] == "plot":
                line = tag[5:]
            elif tag[0:5] == "point":
                point = tag[6:]
        line = self.plots_dict[int(line)]
        point = line.points_dict[int(point)]
        if get_plot:
            return point, line
        else:
            return point

    #################################################################
    # Data handling.                                                #
    #################################################################

    def export_to_csv(self):
        output = StringIO.StringIO()
        for plot in self.plots_list:
            print >> output, plot.label
            for point in plot.points_list:
                if point.additional_data.has_key("export_label"):
                    print >> output, point.additional_data["export_label"], ",", point.xd, "," , point.yd
                else:
                    print >> output, point.xd, "," , point.yd
            print >> output, ""
        contents = output.getvalue()
        output.close()
        return contents

###################################################################################################
# Classes to represent the elements which will appear on the plotting canvas.                     #
###################################################################################################

class Custom_plot:
    def __init__(self, x_data, y_data, plot_id=0, label=None, additional_data=None, color="blue"):
        # Checks the data provided.
        assert(len(x_data) == len(y_data))
        if additional_data:
            assert(len(x_data) == len(additional_data))

        # Initializes the points list.
        self.points_list = [Custom_point(None,None) for i in range(0,len(x_data))]
        self.points_dict = {}

        # Sets the x,y data in the 'data' coordinate system.
        self.set_data_coords(x_data, y_data)
        self.set_additional_data(additional_data)
        self.filtered_points_list = filter(lambda v: v.xd != None and v.yd != None, self.points_list)

        # Sets some attributes of the plot.
        self.id = plot_id
        self.segments_list = []
        self.set_segments()
        self.color = color
        self.label = label


    def get_points(self, exclude_none_values=False):
        if exclude_none_values:
            return self.filtered_points_list
        else:
            return self.points_list

    def set_segments(self):
        segments_points = []
        current_segment = []
        for point in self.points_list:
            if point.yd != None:
                current_segment.append(point)
            else:
                segments_points.append(current_segment)
                current_segment = []
        segments_points.append(current_segment)
        for segment in segments_points:
            if segment != []:
                self.segments_list.append(Custom_plot_segment(segment))


    # Sets data.
    def set_data_coords(self, x_data, y_data):
        self.set_coords(x_data, y_data, system="d")

    def set_canvas_coords(self, x_canvas_coords, y_canvas_coords):
        self.set_coords(x_canvas_coords, y_canvas_coords, system="c")

    def set_coords(self, x_coords, y_coords, system="d"):
        for point,x,y in zip(self.get_points(),x_coords,y_coords):
            setattr(point, "x"+system, x)
            setattr(point, "y"+system, y)

    def set_additional_data(self, additional_data):
        for point, data in zip(self.get_points(), additional_data):
            point.additional_data = data


    # Gets data.
    def get_canvas_coords(self, component="xy", return_mode="lists"):
        return self.get_coords(component=component, system="c", return_mode=return_mode)

    def get_data_coords(self, component="xy", return_mode="lists"):
        return self.get_coords(component=component, system="d", return_mode=return_mode)

    def get_coords(self, component="xy", system="d", return_mode="lists"):
        if component == "xy":
            if return_mode == "lists":
                return map(lambda point: getattr(point, "x"+system), self.get_points()), map(lambda point: getattr(point, "y"+system), self.points_list)
            elif return_mode == "vectors":
                return map(lambda point: (getattr(point, "x"+system), getattr(point, "y"+system)), self.get_points())
        else:
            return map(lambda point: getattr(point, component+system), self.get_points())

    def get_additional_data(self):
        return [point.additional_data for point in self.get_points()]


class Custom_plot_segment:
    def __init__(self, points):
        self.points = points
        self.id = None


class Custom_point:
    def __init__(self, x_data, y_data, color = None):
        self.xd = x_data
        self.yd = y_data

        self.xc = None
        self.yc = None

        self.additional_data = None

        self.color = color
        self.id = None
        self.highlight_id = None


class Custom_line:
    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.id = None


class Custom_label:
    def __init__(self, text, coordinates):
        self.coordinates = coordinates
        self.text = text
        self.id = None


###################################################################################################
# Other functions.                                                                                #
###################################################################################################
def float_range(x_min, x_max, increment):
    if has_numpy:
        a = numpy.arange(x_min, x_max, increment)
        # a = numpy.linspace(x_min, x_max,10,endpoint=False)
        return list(a)
    else:
        range_f = []
        x = x_min
        while x < x_max:
            range_f.append(x)
            x += increment
        return range_f

def custom_linspace(x_min, x_max, steps, endpoint=True):
    if has_numpy:
        a = numpy.linspace(x_min, x_max, steps, endpoint=endpoint)
        return a.tolist()
    else:
        return float_range(x_min, x_max+(x_max-x_min)/(steps-1), (x_max-x_min)/(steps-1))


###################################################################
# Functions needed to build equally spaced intervals on the axes. #
###################################################################
def get_order_of_magnitude(n):
    f = '%e' % n
    if "+" in f:
        e = int(f.split("+")[-1])
    else:
        e = -int(f.split("-")[-1])
    return e

def get_increment(order_of_magnitude):
    if order_of_magnitude == 0:
        return 1
    else:
        return (10**(order_of_magnitude))/2

def round_min(n, order_of_magnitude, to_half = False):
    rn = round(n, -order_of_magnitude)
    if rn > n and to_half:
        rn -= get_increment(order_of_magnitude)
    return rn

def round_max(n, order_of_magnitude, to_half = False):
    rn = round(n, -order_of_magnitude)
    if rn < n and to_half:
        rn += get_increment(order_of_magnitude)
    return rn


def get_plotting_interval(min_val, max_val, mode="use_min_and_max"): # "expand", "use_min_and_max"
    """
    Gets a min and a max value, and returns:
        - a list of equally spaced numbers comprised between the too limits.
        - a list of strings to be used as labels for the numbers in the above list.
    """
    # If 'True', some information regarding the linspace generation will be printed.
    out = False
    # Gets the order of magnitude of the difference between min and max.
    d = max_val - min_val
    dmao = get_order_of_magnitude(d)
    # ---
    if out:
        print "#########################"
        print "# Defining a new plotting interval."
        print "Initial rounding of:",  min_val, max_val
        print "Difference:", d
        print "Difference order of magnitude:", dmao
    # ---

    #-----------------------------------------------------------------------------------------------
    # In this mode the limits will be 'expanded', that is, the min and max will be respectively    -
    # decreased and increased until a nice human-readable linspace can be generated. The results   -
    # are supposed to look like the ones matplotlib generates when initially displaying a plot in  -
    # an without any user supplied parameters (such as axes limits, thicks, etc...).               -
    #-----------------------------------------------------------------------------------------------
    if mode == "expand":
        # ---
        if out:
            print "Expanding the initial limits to find an equally spaced interval."
        # ---
        cmi = round_min(min_val, dmao, to_half=True)
        cma = round_max(max_val, dmao, to_half=True)
        # Transforms the numbers in integers if the order of magnitude their difference is less than 0.
        if dmao < 0:
            cma = cma * 10**(-dmao)
            cmi = cmi * 10**(-dmao)
            tdmao = 0
            cmi = round_min(cmi, tdmao, to_half=True)
            cma = round_max(cma, tdmao, to_half=True)
            # ---
            if out:
                print "Adjusting rounding for oom < 0:", cmi,cma, tdmao
            # ---
        else:
            tdmao = dmao

        # Finds the intervals.
        steps = None
        delta = None
        tested_original = False
        decrease_min_value = True
        counter = 0

        while not steps:
            # Try to find a number of intervals bewteen the min and max value in which the limits are
            # non float numbers.
            for i in (4,5,6,7,8,9,10):
                r = (cma-cmi)%i
                # ---
                if out:
                    print "Try:", r, cma, cmi, (cma-cmi)%i
                # ---
                if r == 0:
                    # ---
                    if out:
                        print "Solution found at cycle:", counter
                    # ---
                    steps = i
                    delta = (cma-cmi)/i
                    break
            if steps:
                break
            # If the solution was not found in the first cycle, then start to change the max and min
            # values of the axes, until a solution is found.
            if decrease_min_value:
                if min_val > 0 and cmi == 0:
                    pass
                else:
                    cmi -= 10**tdmao
                decrease_min_value = False
            else:
                cma += 10**tdmao
                decrease_min_value = True
            counter += 1
            # Breaks the search at cycle 1000: something must have gone wrong.
            if counter > 1000:
                break

        # Reconvert the min and max values in floats, if the order of magnitude of their difference was
        # less than 0.
        if dmao < 0:
            cma = cma / 10**(-dmao)
            cmi = cmi / 10**(-dmao)
            delta = delta / 10**(-dmao)

        # Builds a linspace using the min and max values found above.
        numbers_linspace = custom_linspace(cmi, cma, steps+1)

        # ---
        if out:
            print "#########################"
            print "Results, obtained in cycles:", counter
            print "Adjusted Min and Max:", cmi, cma
            print "Steps:", steps, "Delta:", delta
            print "Linspace:", numbers_linspace
        # ---

    #-----------------------------------------------------------------------------------------------
    # In this mode the min and max values will be mantained and an equally spaced interval         -
    # overlapping the min and max limits will be found.                                            -
    #-----------------------------------------------------------------------------------------------
    elif mode == "use_min_and_max":
        # ---
        if out:
            print "Finding an equally spaced interval overlapping the min and max limits."
        # ---

        # Round the min and max values.
        cmi = round_max(min_val, dmao)
        cma = round_min(max_val, dmao)

        # Generate an initial linspace. Since the min and max have been rounded, some values of the
        # linspace might be found outside of the original limits.
        c = 0
        cc = cmi
        numbers_linspace = []
        while cc < cma:
            numbers_linspace.append(cc)
            cc += 10**dmao
            c += 1
            if c > 100:
                break
        numbers_linspace.append(cma)

        # ---
        if out:
            print "Rounded Min and Max:", cmi, cma, cmi == cma
            print "Initial linspace:", numbers_linspace
        # ---

        if len(numbers_linspace) > 1:
            # ---
            if out:
                print "Adjusting to obtain at least 4 thicks."
            # ---
            filtered_numbers_linspace = filter(lambda v: v >= min_val and v <= max_val, numbers_linspace)
            filtered_numbers_set = set(filtered_numbers_linspace)
            c = 0
            # Try to produce a linspace with at least 4 equally spaced numbers.
            while len(filtered_numbers_set) < 4:
                # Add extra numbers in the middle of the elements of the current linspace.
                for l0, l1 in zip(numbers_linspace[0:], numbers_linspace[1:]):
                    n = l0+(l1-l0)/2.0
                    numbers_linspace.insert(numbers_linspace.index(l0)+1, n)
                    filtered_numbers_set.add(n)
                delta = numbers_linspace[1] - numbers_linspace[0]
                # Add extra thicks if necessary, in order to fill the space between the min and max
                # values of the linspace just generated and the min and max limits of the interval.
                if delta != 0:
                    c1 = 0
                    while min(filtered_numbers_set)-delta >= min_val:
                        filtered_numbers_set.add(min(filtered_numbers_set)-delta)
                        c1 += 1
                        if c1 > 500:
                            break
                    c1 = 0
                    while max(filtered_numbers_set)+delta <= max_val:
                        filtered_numbers_set.add(max(filtered_numbers_set)+delta)
                        c1 += 1
                        if c1 > 500:
                            break
                # When dealing with numbers too small (10^-18), the cycle will have to be stopped,
                # beacause there will be problems finding building the linspace.
                c += 1
                if c > 5:
                    break
            filtered_numbers_linspace = sorted(list(filtered_numbers_set))
            numbers_linspace = filtered_numbers_linspace
            numbers_linspace = filter(lambda v: v >= min_val and v <= max_val, numbers_linspace)

    # Rounds numbers in the linspace which have an order of magnitude smaller than order of
    # magnitude of the difference between the original min and max values.
    rounded_numbers_linspace = []
    for n in numbers_linspace:
        if get_order_of_magnitude(n) < dmao:
            f = "%0."+str(abs(dmao))+"f"
            rounded_numbers_linspace.append(float(f % n)) # numbers_linspace[i] = float(f % n)
        else:
            rounded_numbers_linspace.append(n)
    # ---
    if out:
        print "Formatted linspace",rounded_numbers_linspace
        print "Delta:", numbers_linspace[1] - numbers_linspace[0]
    # ---

    # This will be populated with labels for each number of the linspace.
    strings_linspace = []
    # Maximum number of significant figures allowed in labels.
    msf = 3
    for r in numbers_linspace:# rounded_numbers_linspace:
        rf = None
        # Removes the trailing decimal zeroes from integers.
        if r%1 == 0:
            rf = str(r).split(".")[0]
        else:
            ooms = set([get_order_of_magnitude(i) for i in numbers_linspace])
            if len(ooms) > 1:
                if max(ooms) - min(ooms) > 5 and get_order_of_magnitude(r) == min(ooms):
                    r = round(r,-max(ooms))
            rf = str(r)
        # Converts a '0.0' lable in just '0'.
        if r == 0:
            rf = "0"
        strings_linspace.append(rf)

    if mode == "use_min_and_max" and min_val != cmi:
        numbers_linspace.insert(0,min_val)
        strings_linspace.insert(0,"")
    if mode == "use_min_and_max" and max_val != cma:
        numbers_linspace.append(max_val)
        strings_linspace.append("")

    return numbers_linspace, strings_linspace


###################################################################################################
# Builds a plot showing a distance tree.                                                          #
###################################################################################################

# The following code has been adapted from the 'draw' method of the 'Phylo' module of Biopython.
# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
'''
                 Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
'''
def draw_tree(tree, plotting_window_parent=None, label_func=str, do_show=True, show_confidence=True):
    """
    Draws a distance tree contained in a 'Tree' class object of Bio.Phylo using PyMod plotting
    engine. It mainly uses the algorithm implemented in the 'draw' method of 'Phylo'.
    """
    #############################################################################
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Layout
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights
    #############################################################################

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)

    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    xlim = (-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ylim = (0.2, max(y_posns.values()) + 0.8) # (max(y_posns.values()) + 0.8, 0.2)

    # Aesthetics
    plot_title = "Distance tree"
    if hasattr(tree, 'name') and tree.name:
        plot_title = tree.name

    # Builds the plot window and initializes the plotting area.
    cp = Custom_plot_window(plotting_window_parent, title=plot_title)
    cp.build_plotting_area(use_plotting_controls=False, x_label_text="Branch length", y_label_text="Taxum", hide_y_thicks=True)
    cp.show(x_limits=xlim, y_limits=ylim)

    #############################################################################
    def draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0,
                         color='black', lw='.1'):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == 'horizontal':
            horizontal_linecollections.append([x_start, y_here, x_here, y_here])

        elif orientation == 'vertical':
            vertical_linecollections.append([x_here, y_bot, x_here, y_top])


    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, 'color') and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, 'width') and clade.width is not None:
            lw = 1 # GX: clade.width * plt.rcParams['lines.linewidth']
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw)

        #----------------------------------------------------
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            cp.draw_label(' %s' % label, [x_here, y_here])
        # Add label above the branch (optional)
        show_confidence = True
        if show_confidence:
            def format_branch_label(clade):
                if hasattr(clade, 'confidences'):
                    # phyloXML supports multiple confidences
                    return '/'.join(conf2str(cnf.value)
                                    for cnf in clade.confidences)
                if clade.confidence:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
        conf_label = format_branch_label(clade)
        if conf_label:
            cp.draw_label(conf_label, [0.5 * (x_start + x_here), y_here])
        #----------------------------------------------------

        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)
    #############################################################################

    draw_clade(tree.root, 0, 'k', 1)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for coords in horizontal_linecollections:
        print coords
        cp.draw_line(coords)
    for coords in vertical_linecollections:
        cp.draw_line(coords)


###################################################################################################
# Shows a simple test.                                                                            #
###################################################################################################

if __name__ == "__main__":
    def prova(point, plot):
        print "Test: ", point.xd, point.yd, point.additional_data, plot.label
    root = Tk()
    cp = Custom_plot_window(root, title="A test plotting window.")
    cp.build_plotting_area(message_bar_initial_text="Click on points to highlight residues in PyMOL: ", y_label_text="DOPE profile", on_click_action=prova)
    # vals = list(numpy.random.uniform(-0.5,0.5,200)*100)
    n = 5
    for p in range(0,n):
        l = numpy.random.randint(5, 1200) # 250
        m =  numpy.random.rand() * numpy.random.randint(-100, 100)*0.001 # -0.035
        std =  numpy.random.rand() * numpy.random.randint(-100, 100)*0.001 # 0.01
        vals = std * numpy.random.randn(l) + m # 0.0001, 0.01
        vals = list(vals) # vals = range(0,100)
        for i in numpy.random.randint(low=0,high=l,size=10):
            for j in range(5):
                vals.insert(i,None)
        cp.add_plot(None,list(vals), label="plot %s" % (numpy.random.rand(1)[0]))
    cp.show(x_limits="use_data_limits", y_limits="use_data_limits")

    root.mainloop()
