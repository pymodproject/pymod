from Tkinter import *
from tkFileDialog import *
import tkMessageBox
import tkFont
import Pmw
import os
import sys

# Provides functionalities to some widgets.
import pymod_sequence_manipulation as pmsm

###################################################################################################
# GENERAL STYLES.                                                                                 #
###################################################################################################

# Defines the fixed width font to use in the plugin.
fixed_width_font = None
if sys.platform == "darwin":
    fixed_width_font = "monospace"
else:
    fixed_width_font = "courier"

widgets_background_color = "black"
inactive_entry_bg = "gray"

# Labels.
label_style_0 = {"font": "comic 12", # "height": 2,
                 "background" : widgets_background_color,
                 "fg": "white", "pady": 2}

label_style_1 = {"font": "comic 12",
                 # "height" : 1,
                 "background": widgets_background_color,
                 "fg": "red", "padx": 8}

label_style_2 = {"font": "comic 10", "height": 1,
                 "background": widgets_background_color,
                 "fg": "white", "padx": 10}

small_label_style = {"background": widgets_background_color,
                     "fg":'white', "anchor":"nw",
                     "justify":"left"}

# Buttons.
button_style_1 = {"relief": "raised","borderwidth":"3", "bg":"black", "fg":"white", "pady" : 5} # Submit button.
button_style_2 = {"relief": "raised","borderwidth":"3", "bg":"black", "fg":"white", "pady" : 0}
avdanced_button_style = {"borderwidth": 1, "bg": "black", "fg": "white"}

# Radiobuttons.
radiobutton_style_1 = {"foreground" : "white", "anchor" : "center",
                       "indicatoron":0, "highlightbackground":'black',
                       "selectcolor" : "red", "bg" : 'darkgrey',"padx" : 3,"pady" : 3}
# Pack options.
pack_options_1 = {"side": "top", "padx": 10, "pady": 10, "anchor": "w"}

# ---
# Modeling window styles.
# ---
modeling_window_title_style = {"font": "comic 12", "height": 1,
                        "background":widgets_background_color,
                        "fg":'red', "borderwidth" : 0,
                        "padx" : 20, "pady" : 7}

template_title_options = {"font": tkFont.Font(family="comic",size=10,slant="italic"), "background": widgets_background_color, "fg":'red', "anchor":"w"}

modeling_options_sections_style = {"font": "comic 11", "background": widgets_background_color,
        "fg":'white', # orange, wheat1, orange red
        "anchor":"w", "padx" : 30, "pady" : 7}

modeling_window_option_style = {"font": "comic 10",
                                "background": widgets_background_color,
                                "fg":'red', "anchor":"w"}

modeling_window_explanation = {"font": "comic 8",
                               "background": widgets_background_color,
                               "fg":'white', "anchor":"nw",
                               "padx": 30, "justify":"left"}
modeling_window_explanation_padless = {"font": "comic 8",
                               "background": widgets_background_color,
                               "fg":'white', "anchor":"nw",
                               "justify":"left"}

modeling_window_rb_big = {"bg":widgets_background_color,
                          "highlightbackground":widgets_background_color,
                          "fg":"red", "font":"comic 10"}

modeling_window_rb_small = {"bg":widgets_background_color,
                            "fg":"white", "selectcolor": "red",
                            "highlightbackground":'black'}

# Checkbuttons.
modeling_window_checkbutton = {"background": widgets_background_color,
                               "foreground": "white", "selectcolor": "red",
                               "highlightbackground": 'black'}

# Frames.
target_box_style = {"background":'black', "bd":1, "relief":GROOVE, "padx":10, "pady":10}


###################################################################################################
# CLASSES FOR WIDGETS USED THROUGHOUT THE PLUGIN MULTIPLE TIMES.                                  #
###################################################################################################

class PyMod_gui_style_mixin:
    pass


class PyMod_frame(Frame):
    """
    A class for frames created in the PyMod GUI.
    """
    def __init__(self, parent = None, **configs):
        Frame.__init__(self, parent, background = widgets_background_color, **configs)


class PyMod_base_window(Toplevel):
    """
    A class for a base window created in PyMod.
    """

    def __init__(self, parent = None, title = "PyMod window", freeze_parent=True, **configs):

        Toplevel.__init__(self, parent, **configs)

        self.title("<<" + title + ">>")

        # Frame that occupies all the window.
        self.main_frame = PyMod_frame(self)
        self.main_frame.pack(expand = YES, fill = BOTH)

        if freeze_parent:
            try:
                self.grab_set()
            except:
                pass


class PyMod_tool_window(PyMod_base_window):
    """
    A class to build "Tool Windows", windows that are used in PyMod to launch specific algorithms.
    These windows have three main frames, an upper one (containing a description of what the window
    is used for), a middle one (containing the options needed to set the parameters of some
    algorithm) and a lower on (containing a "SUBMIT" button used to launch the algorithm).
    """
    def __init__(self, parent = None, title = "PyMod window", upper_frame_title="Here you can...", submit_command=None, with_frame=False, pack_options=None , **configs):

        PyMod_base_window.__init__(self, parent, title, **configs)

        # Builds the upper frame with the title.
        self.upperframe = PyMod_frame(self.main_frame, borderwidth=5, relief='groove', pady=15)
        self.upperframe.pack(side = TOP, expand = NO, fill = X, ipadx = 3, ipady = 3, pady=15)
        self.upperframe_title=Label(self.upperframe, text = upper_frame_title, **label_style_0)
        self.upperframe_title.pack(fill="x")

        # Builds the middle frame where the tool's options are going to be displayed.
        if with_frame:
            self.mid_scrolled_frame = Pmw.ScrolledFrame(self.main_frame,
                hull_bg='black', frame_bg='black',
                usehullsize = 0, borderframe = 0, vscrollmode='dynamic', hscrollmode='none',
                hull_borderwidth = 0, clipper_bg='black')
            self.mid_scrolled_frame.pack(side= TOP, fill = BOTH, anchor="w", expand = 1)
            # This is the actual Frame where the content of the tab is going to be packed.
            self.midframe = self.mid_scrolled_frame.interior()
            self.midframe.configure(padx = 5, pady = 5)
        else:
            self.midframe = PyMod_frame(self.main_frame)
            self.midframe.pack(side = TOP, fill = BOTH, anchor="w", ipadx = 5, ipady = 5)

        # Builds the lower frame of the modeling window where the "SUBMIT" button is.
        self.lowerframe = PyMod_frame(self.main_frame)
        self.lowerframe.pack(side = BOTTOM, expand = NO, fill = Y, anchor="center", ipadx = 5, ipady = 5)
        self.submit_button=Button(self.lowerframe, text="SUBMIT", command=submit_command, **button_style_1)
        self.submit_button.pack(pady=10)

        # Define the way in which the window widgets are going to be packed.
        self.pack_options = None
        if pack_options != None:
            self.pack_options = pack_options
        else:
            self.pack_options = pack_options_1

        # A list which is going to contain the widgets that will be aligned using the
        # 'Pmw.alignlabels()' and 'align_input_widgets_components()' functions.
        self.widgets_to_align = []

        # A list of widgets which is going to be hidden until users press the "Advanced" button in
        # oreder to display the window's advanced options.
        self.advanced_widgets = []
        self.showing_advanced_widgets = False

        # A list of widgets whose input is going to be validated once the 'SUBMIT' button is pressed.
        self.widgets_to_validate = []


    def add_widget_to_align(self, widget):
        self.widgets_to_align.append(widget)

    def align_widgets(self, input_widget_width=10):
        Pmw.alignlabels(self.widgets_to_align, sticky="nw")
        align_input_widgets_components(self.widgets_to_align, input_widget_width)


    def add_advanced_widget(self, widget):
        self.advanced_widgets.append(widget)

    def show_advanced_button(self):
        # Pads a little the button towards the window's center.
        button_pack_options = self.pack_options.copy()
        if button_pack_options.has_key("padx"):
            button_pack_options["padx"] += 5
        self.advance_options_button = Button(self.midframe, text="Show Advanced Options", command=self.toggle_advanced_options,**avdanced_button_style)
        self.advance_options_button.pack(**button_pack_options)

    def toggle_advanced_options(self):
        if self.showing_advanced_widgets:
            for w in self.advanced_widgets:
                w.pack_forget()
            self.showing_advanced_widgets = False
            self.advance_options_button.configure(text="Show Advanced Options")
        else:
            for w in self.advanced_widgets:
                w.pack(**self.pack_options)
            self.showing_advanced_widgets = True
            self.advance_options_button
            self.advance_options_button.configure(text="Hide Advanced Options")


    def add_widget_to_validate(self, widget):
        self.widgets_to_validate.append(widget)

    def get_widgets_to_validate(self):
        widgets = []
        for w in self.widgets_to_validate:
            if not w in self.advanced_widgets:
                widgets.append(w)
            else:
                if self.showing_advanced_widgets:
                    widgets.append(w)
        return widgets


class PyMod_radioselect(Pmw.RadioSelect):
    """
    Class for custom Pmw.RadioSelect widgets.
    """
    def __init__(self, parent = None, label_style=None, **configs):
        Pmw.RadioSelect.__init__(self, parent,buttontype = 'radiobutton',
                                 orient = 'vertical', labelpos = 'wn',
                                 pady=0, padx=0,
                                 labelmargin=5,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("frame").configure(background=widgets_background_color)
        self.component("hull").configure(background=widgets_background_color)

        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)
        self.buttons_list = []


    def add(self, componentName, **configs):
        """
        Override the "add" method in order to change the Buttons appearance according to PyMod
        style.
        """
        widget = Pmw.RadioSelect.add(self, componentName, **configs)
        widget.configure(**radiobutton_style_1)
        self.buttons_list.append(widget)

    def set_input_widget_width(self, widget_width):
        for b in self.buttons_list:
            b.configure(width = widget_width)


class PyMod_entryfield(Pmw.EntryField):
    """
    Class for custom Pmw.EntryField widgets.
    """
    def __init__(self, parent = None, label_style=None, run_after_selection=None, **configs):
        Pmw.EntryField.__init__(self, parent,
                                 labelpos = 'wn',
                                 labelmargin=5,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("hull").configure(background=widgets_background_color)
        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)
        self.run_after_selection = run_after_selection

    def set_input_widget_width(self, widget_width):
        self.component("entry").configure(width = widget_width)


class PyMod_path_entryfield(PyMod_entryfield):
    """
    Class for a custom entryfield accompanied by a button to choose a path on the user's machine.
    """
    def __init__(self, parent = None, label_style=None, path_type="file", askpath_title = "Search for a path", file_types="", **configs):
        PyMod_entryfield.__init__(self, parent, label_style, **configs)

        self.interior = self.component('hull') # self.interior()
        self.path_type = path_type
        self.file_types = file_types
        self.askpath_title = askpath_title

        self.choose_path_button = Button(self.interior, text="Browse", command=self.choose_path, **button_style_2)
        self.choose_path_button.grid(column=3,row=2, padx=(15,0))

        self.component("entry").configure(readonlybackground=inactive_entry_bg)


    def choose_path(self):
        """
        Called when users press the 'Browse' button in order to choose a path on their system.
        """
        current_path = self.getvalue()
        new_path = None

        # Lets users choose a new path.
        if self.path_type == "file":
            new_path = askopenfilename(title = self.askpath_title,
                initialdir=os.path.dirname(current_path),
                initialfile=os.path.basename(current_path), parent = get_parent_window(self), filetypes = self.file_types)

        elif self.path_type == "directory":
            new_path = askdirectory(title = self.askpath_title, initialdir=os.path.dirname(current_path), mustexist = True, parent = get_parent_window(self))

        # Updates the text in the Entry with the new path name.
        if new_path:
            self.clear()
            self.setvalue(new_path)

        if hasattr(self.run_after_selection, "__call__"):
            self.run_after_selection()


class PyMod_combobox(Pmw.ComboBox):
    """
    Class for custom combobox widgets.
    """
    def __init__(self, parent = None, label_style=None, **configs):
        Pmw.ComboBox.__init__(self, parent,
                                 labelpos = 'wn',
                                 labelmargin=5,history = 0,
                                 **configs)
        # Configure the widgets component to change their appearance according to PyMod style.
        self.component("hull").configure(background=widgets_background_color)
        if label_style == None:
            label_style = label_style_1
        self.component("label").configure(**label_style)

        # Configure other properties of the combobox.
        self.component("entryfield").component("entry").configure(state='readonly',
            readonlybackground= "white", fg="black", bg="white")
        self.selectitem(0)

    def set_input_widget_width(self, widget_width):
        self.component("entryfield").component("entry").configure(width = widget_width)


class PyMod_dialog(Pmw.MessageDialog):
    def __init__(self, parent, **configs):
        Pmw.MessageDialog.__init__(self, parent, command = self.dialog_state,**configs)
        self.wait_state = True
        self.val = None

    def dialog_state(self,val):
        self.withdraw()
        self.val = val
        self.wait_state = False
        return self.val

    def get_dialog_value(self):
        if self.wait_state:
            self.after(100, self.get_dialog_value)
        val = self.val
        # Resets the value to default.
        self.wait_state = True
        self.val = None
        return val


###################################################################################################
# CLASSES FOR WIDGETS USED IN SPECIFIC PARTS OF THE PyMod GUI.                                    #
###################################################################################################

#####################################################################
# MODELLER options in the PyMod options window.                     #
#####################################################################

class Use_importable_modeller_radioselect(PyMod_entryfield):

    def __init__(self, parent = None, label_style=None, importable_modeller=False, initial_value=None, **configs):
        PyMod_entryfield.__init__(self, parent, label_style, **configs)
        self.interior = self.component('hull') # self.interior()
        self.component("entry").configure(state='readonly', readonlybackground= "black", bd=0, relief=FLAT, fg="white", bg="black")
        self.importable_modeller = importable_modeller
        # If MODELLER libs are importable, then build some radiobuttons to choose whether to us it.
        self.gui_var = StringVar()
        self.gui_var.set(str(initial_value))
        if self.importable_modeller:
            rbs = radiobutton_style_1.copy()
            rbs["padx"] = 0
            use_radiobutton = Radiobutton(self.interior, text="Use", variable=self.gui_var, value="True", width=8, command=self.run_after_selection, **rbs)
            use_radiobutton.grid(row=2, column=3, sticky = "w",padx=(15,0), pady=(0,3))
            dont_use_radiobutton = Radiobutton(self.interior, text="Don't use", variable=self.gui_var, value="False", width=8, command=self.run_after_selection, **rbs)
            dont_use_radiobutton.grid(row=3, column=3, sticky = "w",padx=(15,0))

    def getvalue(self):
        return self.gui_var.get()


class Modeller_exec_entryfield(PyMod_path_entryfield):

    def __init__(self, parent = None, label_style=None, path_type="file", askpath_title = "Search for a path", file_types="", **configs):
        PyMod_path_entryfield.__init__(self, parent, label_style, path_type=path_type, askpath_title=askpath_title, file_types=file_types, **configs)
        self.not_necessary_label = Label(self.interior, text="Not necessary",**small_label_style)

    def show_path_selector(self, path_to_show):
        self.component("entry").configure(state=NORMAL)
        self.not_necessary_label.grid_forget()
        self.choose_path_button.grid(column=3,row=2, padx=(15,0))

    def hide_path_selector(self):
        self.component("entry").configure(state=NORMAL)
        self.choose_path_button.grid_forget()
        self.not_necessary_label.grid(column=3,row=2, padx=(15,0))
        self.component("entry").configure(state="readonly")


#####################################################################
# Modeling window classes.                                          #
#####################################################################

class Structure_frame:
    """
    A class to construct the template selection frame and to store all their tkinter widgets and
    information.
    """
    labels_width = 14
    template_options_style = modeling_window_option_style.copy()
    template_options_style.update({"width": labels_width})
    frames_padding = 7

    def __init__(self, pymod_object, structure_pymod_element,target_pymod_element,target_widget,structure_number, modeling_cluster_number):
        # These will contain a Structure type object.
        self.structure_pymod_element = structure_pymod_element
        self.target_pymod_element = target_pymod_element
        # The widget in which to grid the structure Frame.
        self.target_widget = target_widget
        # The int value that is passed in the for cycle in which the Structure_frame objects are
        # constructed. Identifies different Structure_frame objects.
        self.id = structure_number
        self.mc_id = modeling_cluster_number # This is the id of the modeling cluster containing a structure frame.
        # This is needed to check what is the state of radiobutton for using hetres. If it is on,
        # then this value should be 1 (by default it is 1 because of the default state of the
        # radiobutton), when it is off, this vaule should be 0.
        self.hetres_radiobutton_state = 1
        self.pymod_object = pymod_object


    def build_frame(self):
        """
        Builds a frame for each template structure and all its options.
        """
        self.structure_frame = Frame(self.target_widget, **target_box_style)
        self.structure_frame.pack(anchor="w",padx=30,pady=(0,5))
        self.build_use_structure_frame()
        # self.build_sequence_limit_frame()
        self.build_hetres_frame()
        self.build_water_frame()


    def build_use_structure_frame(self):
        """
        Builds a Frame which will contain the the checkbox for using the structure as a template.
        """
        # Use-structure frame
        self.use_structure_frame = Frame(self.structure_frame, background='black', pady = Structure_frame.frames_padding)
        self.use_structure_frame.grid(row=0, column=0,sticky = "w")

        # Label for the structure
        self.template_title_lab = Label(self.use_structure_frame, text= "", **template_title_options)
        self.template_title_lab.pack(side = TOP, anchor="w",pady = (0, Structure_frame.frames_padding))

        self.lab=Label(self.use_structure_frame, text= "Use as Template: ",**Structure_frame.template_options_style)
        self.lab.pack(side = LEFT)

        # Checkbutton for using the structure as a template.
        self.use_as_template_var = IntVar()
        # Avoids some problems templates with very long names.
        # For really long names it will print something like: sp_P62987...Chain:A
        template_name = self.structure_pymod_element.my_header[0:-8]
        if len(template_name) < 10:
            template_name = self.structure_pymod_element.my_header
        else:
            template_name = self.structure_pymod_element.my_header[0:10]+"..."+self.structure_pymod_element.my_header[-7:]
        # Shows the identity % between the two aligned sequences.
        identity = pmsm.compute_sequence_identity(self.target_pymod_element.my_sequence, self.structure_pymod_element.my_sequence)
        checkbox_text = template_name + " (id: " + str(identity) + "%)"
        self.chk = Checkbutton(self.use_structure_frame, text=checkbox_text, variable=self.use_as_template_var,command = self.click_on_structure_checkbutton, **modeling_window_checkbutton)
        self.chk.pack(side = LEFT)

        self.template_title_lab["text"] = "Options for template: " + template_name


    def click_on_structure_checkbutton(self):
        """
        This is called when the checkbutton to use the structure as a template is pressed. If users
        want to use hetero-atoms this method will activate the hetres and water checkbuttons, that
        by default are disabled.
        """
        # This is all under the influence of the state of the "use hetatm" radiobutton in the
        # options page.
        if self.hetres_radiobutton_state == 1:
            # The template is "activated", and so also its checkbuttons are.
            if self.use_as_template_var.get() == 1:
                self.activate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.activate_het_checkbuttons()
            # The template is "inactivated", also its checkbuttons are.
            elif self.use_as_template_var.get() == 0:
                self.inactivate_water_checkbutton()
                if self.number_of_hetres > 0:
                    self.inactivate_het_checkbuttons()


    # This is launched when the hetres radiobutton state is changed to "NO".
    def inactivate_het_checkbuttons(self):
        self.use_all_hetres.configure(state=DISABLED)
        self.select_single_hetres.configure(state=DISABLED)
        self.do_not_use_hetres.configure(state=DISABLED)
        for c in self.structure_hetres_checkbuttons:
            c.configure(state=DISABLED)


    # This is launched when the hetres radiobutton state is changed to "YES".
    def activate_het_checkbuttons(self):
        if self.use_as_template_var.get() == 1:
            self.use_all_hetres.configure(state=NORMAL)
            self.select_single_hetres.configure(state=NORMAL)
            self.do_not_use_hetres.configure(state=NORMAL)
            for c in self.structure_hetres_checkbuttons:
                c.configure(state=NORMAL)


    def activate_water_checkbutton(self):
        if self.structure_pymod_element.structure.has_waters():
            if self.use_as_template_var.get() == 1:
                self.water_checkbox.configure(state=NORMAL)


    def inactivate_water_checkbutton(self):
        if self.structure_pymod_element.structure.has_waters():
            self.water_checkbox.configure(state=DISABLED)


    def build_sequence_limit_frame(self):
        """
        Frame for the sequence limits.
        """
        # From-to frame
        self.limits_frame = Frame(self.structure_frame, background='black', pady = Structure_frame.frames_padding)
        self.limits_frame.grid(row=1, column=0,sticky = "w")

        # From label. The width is relative to the font size
        self.from_enf = PyMod_entryfield(self.limits_frame, label_text = "From: ", value = 1,
                                       validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000},
                                       label_style =modeling_window_option_style)
        self.from_enf.component("entry").configure(width = 5)
        self.from_enf.pack(side="left", padx=0)
        # To label
        self.to_enf = PyMod_entryfield(self.limits_frame, label_text = "To: ", value = 1000,
                                     validate = {'validator' : 'integer', 'min' : 1, 'max' : 5000},
                                     label_style= modeling_window_option_style)
        self.to_enf.component("entry").configure(width = 5)
        self.to_enf.pack(side="left", padx=(20,0))


    def build_hetres_frame(self):
        """
        Builds a frame for the Hetero-residues selection.
        """
        # This is going to contain the checkbox states of the HETRES of the structure.
        self.structure_hetres_states  = []
        self.structure_hetres_checkbuttons = []
        # Hetero-residues frame
        self.hetres_frame = Frame(self.structure_frame, background='black', pady = Structure_frame.frames_padding)
        self.hetres_frame.grid(row=2, column=0,sticky = "w")
        # Label
        self.hetres_label = Label(self.hetres_frame, text= "Hetero Residues: ", **Structure_frame.template_options_style)
        self.hetres_label.grid(row=0, column=0, sticky = "nw")
        # Variable for the radiobuttons.
        self.hetres_options_var = IntVar()
        # Counts the hetres of this chain.
        self.number_of_hetres = len(self.structure_pymod_element.structure.hetero_residues)

        if self.number_of_hetres > 0:
            # Radiobuttons for hetres options and their frame.
            self.hetres_options_frame = Frame(self.hetres_frame, background='black')
            self.hetres_options_frame.grid(row=0, column=1, sticky = "nw")
            self.hetres_options_var.set(1)
            self.use_all_hetres_text = "Use all heteroatomic residues (%s)" % (self.number_of_hetres)
            self.use_all_hetres = Radiobutton(self.hetres_options_frame, text=self.use_all_hetres_text, variable=self.hetres_options_var, value=1,background='black', foreground = "white", selectcolor = "red", highlightbackground='black',command=self.hide_select_single_hetres_frame, state=DISABLED) #,command=self.activate_template_dsb_frame)
            self.use_all_hetres.grid(row=0, column=0, sticky = "w")
            # Select single hetres manually.
            self.select_single_hetres = Radiobutton(self.hetres_options_frame, text="Select single heteroatomic residues", variable=self.hetres_options_var, value=2, command=self.show_select_single_hetres_frame, state=DISABLED,**modeling_window_rb_small)
            self.select_single_hetres.grid(row=1, column=0, sticky = "w")
            self.select_single_hetres_frame = Frame(self.hetres_options_frame, background='black')

            # This is needed to count the "rows" used to grid HETRES checkboxes.
            self.hetres_counter = 0
            for hetres in self.structure_pymod_element.structure.hetero_residues:
                # Checkbox for each HETRES.
                single_hetres_state = IntVar()
                # Complete it with the full name.
                checkbox_text = "%s (%s) %s" % (hetres.three_letter_code,hetres.hetres_type, hetres.pdb_position)
                hetres_checkbutton = Checkbutton(self.select_single_hetres_frame, text=checkbox_text, variable=single_hetres_state, state=DISABLED, **modeling_window_checkbutton)
                hetres_checkbutton.grid(row=self.hetres_counter, column=0, sticky = "w",padx=(15,0))
                self.hetres_counter += 1
                # Adds the single HETRES state to a list that contains the ones of the structure.
                self.structure_hetres_states.append(single_hetres_state)
                self.structure_hetres_checkbuttons.append(hetres_checkbutton)

            self.do_not_use_hetres = Radiobutton(self.hetres_options_frame, text="Do not use any heteroatomic residue", variable=self.hetres_options_var, value=3,command=self.hide_select_single_hetres_frame, state=DISABLED, **modeling_window_rb_small)
            self.do_not_use_hetres.grid(row=3, column=0, sticky = "w")

        else:
            self.no_hetres_label = Label(self.hetres_frame, text="No heteroatomic residue found",background='black', foreground = "gray45")
            self.no_hetres_label.grid(row=0, column=1, sticky = "w")
            self.hetres_options_var.set(3)


    def show_select_single_hetres_frame(self):
        self.select_single_hetres_frame.grid(row=2, column=0, sticky = "w")
        self.pymod_object.templates_frame.reposition()


    def hide_select_single_hetres_frame(self):
        self.select_single_hetres_frame.grid_remove()
        self.pymod_object.templates_frame.reposition()


    def build_water_frame(self):
        """
        Builds a frame for letting the user choose to include water molecules in the model.
        """
        # Frame for water
        self.water_frame = Frame(self.structure_frame, background='black', pady = Structure_frame.frames_padding)
        self.water_frame.grid(row=3, column=0,sticky = "w")
        # Label for water
        self.water_label = Label(self.water_frame, text= "Include Water: ", **Structure_frame.template_options_style)
        self.water_label.grid(row=0, column=0, sticky = "w")

        # Checkbox for water
        # Variable with the state for including water molecules
        self.water_state = IntVar()
        self.water_state.set(0)
        if self.structure_pymod_element.structure.has_waters():
            n_water = self.structure_pymod_element.structure.water_molecules_count
            self.text_for_water_checkbox = "%s water molecules" % (n_water)
            self.water_checkbox = Checkbutton(self.water_frame, text=self.text_for_water_checkbox, variable=self.water_state, command= lambda x=self.id: self.click_on_water_checkbutton(x),state=DISABLED, **modeling_window_checkbutton)
            self.water_checkbox.grid(row=0, column=1, sticky = "w")
        else:
            self.no_water_label = Label(self.water_frame, text= "This structure has no water molecules", background='black', fg='gray45', anchor ="w")
            self.no_water_label.grid(row=0, column=1, sticky = "w")


    def click_on_water_checkbutton(self,x):
        """
        When a structure water checkbutton is pressed, this method deselects the water checkbutton of
        all the other structures, because only water from one structure can be used to build the
        model.
        """
        # This works but the use of pymod object should be avoided.
        for sf in self.pymod_object.modeling_clusters_list[self.mc_id].structure_frame_list:
            if sf.id != self.id:
                if sf.structure_pymod_element.structure.has_waters():
                    sf.water_checkbox.deselect()


class Disulfides_frame:
    """
    A class to construct disulfide frame in the modeling window and to store all their information.
    """
    dsb_building_mode_label = modeling_window_option_style.copy()
    dsb_building_mode_label.update({"padx": 20, "pady": 7})

    def __init__(self, pymod_object, target_widget):
        # The widget in which to build the frame.
        self.target_widget = target_widget
        self.main_disulfides_frame = Frame(self.target_widget, background='black')
        self.main_disulfides_frame.grid(row=0, column=0, sticky = "nw",pady=(0,5))
        self.pymod_object = pymod_object


    def check_templates_with_dsb(self):
        """
        Checks if there are some templates with disulfide bridges. It returns True if there is at
        least one template with a dsb.
        """
        self.templates_with_dsb = False
        for mc in self.pymod_object.modeling_clusters_list:
            if mc.has_structures_with_disulfides():
                self.templates_with_dsb = True
                break
        return self.templates_with_dsb


    def build_template_dsb_frame(self):
        """
        Builds the top frame, for the templates disulfides.
        """
        # Label for the name of the target.
        self.target_name_label = Label(self.main_disulfides_frame,text="Disulfide options",**modeling_window_title_style)
        self.target_name_label.grid(row=0, column=0, sticky = "nw")

        # The frame for template disulfides.
        self.template_dsb_frame = Frame(self.main_disulfides_frame, background='black')
        self.template_dsb_frame.grid(row=1, column=0, sticky = "nw")

        # Label for the title.
        self.templates_dsb_label = Label(self.template_dsb_frame, text= "Use template disulfides", **Disulfides_frame.dsb_building_mode_label)
        self.templates_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(0,0))

        # If there are some templates with disulfide bridges.
        if self.check_templates_with_dsb():
            # Label for the information about the use of this feature.
            information = "Include disulfide bridges found in the structures in the Templates page."
            self.template_disulfides_information = Label(self.template_dsb_frame, text= information, **modeling_window_explanation)
            self.template_disulfides_information.grid(row=1, column=0, sticky = "w")
            # Radiobuttons and their frame.
            self.use_template_dsb_var = IntVar()
            # Initialize the radiobuttons.
            self.use_template_dsb_var.set(1)
            # Frame.
            self.use_template_dsb_rad_frame = Frame(self.template_dsb_frame)
            self.use_template_dsb_rad_frame.grid(row=2, column=0, sticky = "w")
            # Radiobuttons.
            self.use_template_dsb_rad1 = Radiobutton(self.use_template_dsb_rad_frame, text="Yes", variable=self.use_template_dsb_var, value=1, padx=20,command=self.activate_template_dsb_frame, **modeling_window_rb_big)
            self.use_template_dsb_rad1.grid(row=0, column=0, sticky = "w")
            self.use_template_dsb_rad2 = Radiobutton(self.use_template_dsb_rad_frame, text="No", variable=self.use_template_dsb_var, value=0, padx=20,command=self.inactivate_template_dsb_frame, **modeling_window_rb_big)
            self.use_template_dsb_rad2.grid(row=0, column=1, sticky = "w")
            # Button for displaying the list of disulfide bridges found in the templates.
            self.toggle_template_frame = Frame(self.template_dsb_frame,bg="black")
            self.toggle_template_frame.grid(row = 3, column = 0,sticky = "w",padx = (30,0),pady = (5,0))
            toggle_template_dsb_text = "List of templates' disulfides (white: conserved in target, gray: not conserved):"
            self.toggle_template_dsb_label = Label(self.toggle_template_frame,text=toggle_template_dsb_text,bg="black", fg="white")
            self.toggle_template_dsb_label.grid(row = 0, column = 0,sticky = "w",padx = (0,10))
            self.toggle_template_dsb_button = Button(self.toggle_template_frame,text="Show",command = self.show_template_dsb,**button_style_1)
            self.toggle_template_dsb_button.grid(row = 0, column = 1,sticky = "w")
            self.build_templates_disulfides_frame()

        # If there aren't templates with disulfide bridges.
        else:
            # Label for the information about the use of this feature.
            information = "There isn't any template with disulfide bridges."
            self.template_disulfides_information = Label(self.template_dsb_frame, text = information, **modeling_window_explanation)
            self.template_disulfides_information.grid(row=1, column=0, sticky = "w")


    # Called when the "Yes" radiobutton of the "Use template disulfide" option is pressed.
    def activate_template_dsb_frame(self):
        self.toggle_template_frame.grid(row = 3, column = 0,sticky = "w",padx = (30,0),pady = (5,0))
        self.pymod_object.disulfides_scrolled_frame.reposition()

    # Called when the "Show" button is pressed to show the list of the dsb of the templates.
    def show_template_dsb(self):
        self.template_disulfides_frame.grid(row=4, column=0,sticky = "w",padx = (30,0),pady = (5,0))
        self.toggle_template_dsb_button.configure(text="Hide",command = self.hide_template_dsb)
        self.pymod_object.disulfides_scrolled_frame.reposition()

    # Called when the "No" radiobutton of the "Use template disulfide" option is pressed.
    # This is also called when the "Yes" radiobutton of the "Automatically build disulfides" is
    # pressed.
    def inactivate_template_dsb_frame(self):
        self.toggle_template_frame.grid_remove()
        self.hide_template_dsb()
        self.pymod_object.disulfides_scrolled_frame.reposition()

    # Called when the "Show" button is pressed to hide the list of the dsb of the templates.
    def hide_template_dsb(self):
        self.template_disulfides_frame.grid_remove()
        self.toggle_template_dsb_button.configure(text="Show",command = self.show_template_dsb)
        self.pymod_object.disulfides_scrolled_frame.reposition()

    def build_templates_disulfides_frame(self):
        """
        Builds the frame for displaying disulfide bridges found in the templates.
        """
        # Frame for template disulfides.
        self.template_disulfides_frame = Frame(self.template_dsb_frame, background='black', bd=1, relief = GROOVE, padx = 15, pady = 10)
        # Build a frame for every modeling cluster which have templates with disulfides.
        for mci,mc in enumerate(filter(lambda x:x.has_structures_with_disulfides(),self.pymod_object.modeling_clusters_list)):
            # A counter to iterate through all the template structures.
            frame_for_cluster_templates_dsb = Frame(self.template_disulfides_frame, background='black')
            frame_for_cluster_templates_dsb.grid(row=mci, column=0,sticky = "w", pady=(0,10))
            target_label = Label(frame_for_cluster_templates_dsb, text= "Template dsb for target " + mc.target_name, background='black', fg='red', anchor ="nw",font = "comic 9")
            target_label.grid(row=0, column=0, sticky = "w")

            for ei,element in enumerate(filter(lambda x : x.structure.has_disulfides(), mc.structure_list)):
                disulfides_counter = 0
                # Frame for each structure.
                structure_frame_for_disulfides = Frame(frame_for_cluster_templates_dsb, background='black')
                structure_frame_for_disulfides.grid(row=ei+1, column=0,sticky = "w", pady=(0,10))
                # Label with the name of the structure.
                disulfides_label = Label(structure_frame_for_disulfides, text = element.my_header, background='black', fg='red', width = 14, anchor ="nw",bd = 0, relief = GROOVE,padx = 0)
                disulfides_label.grid(row=0, column=0, sticky = "w")
                # Begins a for cycle that is going to examine all disulfides bridges of the chain.
                for dsb in element.structure.disulfides:
                    # For now, display only intrachain bridges.
                    if dsb.bridge_type == "intrachain":
                        # Check if there are homologous CYS in the target according to the alignment.
                        # Take the target sequence.
                        target = mc.target.my_sequence
                        # CYS 1.
                        cys1_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys1_seq_number)
                        cys1_target_position = pmsm.get_residue_id_in_gapless_sequence(target,cys1_alignment_position) + 1
                        cys1_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target, dsb.cys1_seq_number)
                        cys1_homologue_residue = target[cys1_alignment_position] # The corresponding residue in the target.
                        # CYS 2.
                        cys2_alignment_position = pmsm.get_residue_id_in_aligned_sequence(element.my_sequence, dsb.cys2_seq_number)
                        cys2_target_position = pmsm.get_residue_id_in_gapless_sequence(target,cys2_alignment_position) + 1
                        cys2_is_conserved = pmsm.find_residue_conservation(element.my_sequence, target,dsb.cys2_seq_number)
                        cys2_homologue_residue = target[cys2_alignment_position] # The corresponding residue in the target.
                        # If both CYS that form the disulfide in the template are conserved in the target.
                        if cys1_is_conserved and cys2_is_conserved:
                            # Prints also if the CYS are conserved in the target according to the
                            # alignment.
                            label_text = "Template: C%s - C%s / Target: C%s - C%s" % (dsb.cys1_pdb_number, dsb.cys2_pdb_number, cys1_target_position, cys2_target_position)
                            disulfide_label = Label(structure_frame_for_disulfides, text=label_text, background='black', foreground = "white")

                        else:
                            label_text = "Template: C%s - C%s / Target: %c%s - %c%s" % (dsb.cys1_pdb_number,dsb.cys2_pdb_number, cys1_homologue_residue, cys1_target_position, cys2_homologue_residue, cys2_target_position)
                            disulfide_label = Label(structure_frame_for_disulfides, text=label_text, background='black', foreground = "gray45")
                        disulfide_label.grid(row=disulfides_counter, column=1, sticky = "w")
                        disulfides_counter += 1


    def build_user_defined_dsb_frame(self):
        """
        Builds the bottom frame, for the user-defined disulfides.
        """
        self.user_defined_dsb_frame = Frame(self.main_disulfides_frame, background='black')
        self.user_defined_dsb_frame.grid(row=2, column=0, sticky = "nw")

        self.user_dsb_label = Label(self.user_defined_dsb_frame, text= "Create new disulfides", **Disulfides_frame.dsb_building_mode_label)
        self.user_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(20,0))

        information = "Define custom disulfide bridges to be included in the model. "
        information += ("NOTE: if the S atoms of\n"+
                       "the two cysteines you selected are going to be located more than 2.5A apart in the\n"+
                       "model, MODELLER will not build the bridge." )

        self.user_disulfides_information = Label(self.user_defined_dsb_frame, text = information, **modeling_window_explanation)
        self.user_disulfides_information.grid(row=1, column=0, sticky = "w")

        # Radiobuttons and their frame.
        self.use_user_defined_dsb_var = IntVar()
        self.use_user_defined_dsb_var.set(0)

        # Frame.
        self.use_user_defined_dsb_rad_frame = Frame(self.user_defined_dsb_frame)
        self.use_user_defined_dsb_rad_frame.grid(row=2, column=0, sticky = "w")

        # Radiobuttons.
        self.use_user_defined_dsb_rad1 = Radiobutton(self.use_user_defined_dsb_rad_frame, text="Yes", variable=self.use_user_defined_dsb_var, value=1,padx = 20,command=self.activate_combo_box_frame, **modeling_window_rb_big)
        self.use_user_defined_dsb_rad1.grid(row=0, column=0, sticky = "w")

        self.use_user_defined_dsb_rad2 = Radiobutton(self.use_user_defined_dsb_rad_frame, text="No", variable=self.use_user_defined_dsb_var, value=0, padx = 20,command=self.inactivate_combo_box_frame, **modeling_window_rb_big)
        self.use_user_defined_dsb_rad2.grid(row=0, column=1, sticky = "w")

        # Frame where comboboxes and buttons for user defined disulfides are going to be placed.
        # This is going to be gridded by the "activate_combo_box_frame()" method below.
        self.user_defined_dsb_combo_box_frame = Frame(self.user_defined_dsb_frame, background='black',pady = 5)

        # This will contain a list of User_dsb_selector objects that will store the information
        # about user defined dsb.
        self.user_dsb_selector_list = []


    def activate_combo_box_frame(self):
        self.user_defined_dsb_combo_box_frame.grid(row=3, column=0,sticky = "nw",padx = (30,0))
        self.pymod_object.disulfides_scrolled_frame.reposition()

    def inactivate_combo_box_frame(self):
        self.user_defined_dsb_combo_box_frame.grid_remove()
        self.pymod_object.disulfides_scrolled_frame.reposition()


    def build_auto_dsb_frame(self):
        """
        Builds a frame to display the option to make Modeller automatically create all dsb of the
        model.
        """
        self.auto_dsb_frame = Frame(self.main_disulfides_frame, background='black')
        self.auto_dsb_frame.grid(row=3, column=0, sticky = "nw",pady=(0,25))

        self.auto_dsb_label = Label(self.auto_dsb_frame, text= "Automatically build disulfides", **Disulfides_frame.dsb_building_mode_label)
        self.auto_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(20,0))

        information = ("MODELLER will build a disulfide for every pair of cysteine if they are sufficently close in\n"+
                       "the model. ")
        information += "NOTE: by using this option you will not be able to use the two options above."

        self.auto_disulfides_information = Label(self.auto_dsb_frame, text= information, **modeling_window_explanation)
        self.auto_disulfides_information.grid(row=1, column=0, sticky = "w")

        # Radiobuttons and their frame.
        self.auto_dsb_var = IntVar()
        self.auto_dsb_var.set(0)

        # Frame.
        self.use_auto_dsb_rad_frame = Frame(self.auto_dsb_frame)
        self.use_auto_dsb_rad_frame.grid(row=2, column=0, sticky = "w")

        # Radiobuttons.
        self.auto_dsb_rad1 = Radiobutton(self.use_auto_dsb_rad_frame, text="Yes", variable=self.auto_dsb_var, value=1, padx = 20,command=self.activate_auto_dsb,**modeling_window_rb_big)
        self.auto_dsb_rad1.grid(row=0, column=0, sticky = "w")

        self.auto_dsb_rad2 = Radiobutton(self.use_auto_dsb_rad_frame, text="No", variable=self.auto_dsb_var, value=0,padx = 20, command=self.inactivate_auto_dsb,**modeling_window_rb_big)
        self.auto_dsb_rad2.grid(row=0, column=1, sticky = "w")


    def activate_auto_dsb(self):
        # Inactivates the "use template dsb" radiobuttons and selects the "No" radiobutton.
        if self.templates_with_dsb:
            self.use_template_dsb_rad2.select()
            self.use_template_dsb_rad1.configure(state=DISABLED)
            self.use_template_dsb_rad2.configure(state=DISABLED)
            self.inactivate_template_dsb_frame()

        # Inactivates the "create new dsb" radiobuttons and selects the "No" radiobutton.
        self.use_user_defined_dsb_rad2.select()
        self.use_user_defined_dsb_rad1.configure(state=DISABLED)
        self.use_user_defined_dsb_rad2.configure(state=DISABLED)

        self.user_defined_dsb_combo_box_frame.grid_remove()
        self.pymod_object.disulfides_scrolled_frame.reposition()

    def inactivate_auto_dsb(self):
        # Reactivates the "use template dsb" and the "create new dsb" radiobuttons.
        if self.templates_with_dsb:
            self.use_template_dsb_rad1.configure(state=NORMAL)
            self.use_template_dsb_rad2.configure(state=NORMAL)

        self.use_user_defined_dsb_rad1.configure(state=NORMAL)
        self.use_user_defined_dsb_rad2.configure(state=NORMAL)
        self.pymod_object.disulfides_scrolled_frame.reposition()


    def build_no_dsb_frame(self):
        """
        Builds a frame that is displayed if the target sequence has less than 2 cys.
        """
        self.no_dsb_frame = Frame(self.main_disulfides_frame, background='black')
        self.no_dsb_frame.grid(row=1, column=0, sticky = "nw")

        self.no_dsb_label = Label(self.no_dsb_frame,text= "No disulfide bridge can be built.", **modeling_window_title_style)
        self.no_dsb_label.grid(row=0, column=0, sticky = "nw", pady=(0,0))

        information = "No target sequence has at least two CYS residues needed to form a bridge."
        self.no_disulfides_information = Label(self.no_dsb_frame, text= information, **modeling_window_explanation)
        self.no_disulfides_information.grid(row=1, column=0, sticky = "w")


    def build_modeling_cluster_users_dsb_frame(self, modeling_cluster, modeling_cluster_index):
        """
        Builds a frame where are going to be gridded a series of frames (one for each modeling
        cluster) in order to let the user define additional disulfide bridges for each target.
        """
        modeling_cluster_custom_dsb_frame = Frame(self.user_defined_dsb_combo_box_frame, background='black', bd=1, relief=GROOVE, padx = 10, pady = 10)
        modeling_cluster_custom_dsb_frame.grid(row=modeling_cluster_index, column=0,sticky = "nw",pady = (0,5))
        label_text = ""
        if modeling_cluster.target_with_cys:
            label_text = "Select two CYS for target " + modeling_cluster.target_name
        else:
            label_text = "Target " + modeling_cluster.target_name + " doesn't have at least two CYS residues."
        modeling_cluster_custom_dsb_label = Label(modeling_cluster_custom_dsb_frame,font = "comic 9", text=label_text, bg="black", fg= "red")
        modeling_cluster_custom_dsb_label.grid(row=0, column=0,sticky = "nw")
        uds = User_dsb_selector(self.pymod_object, modeling_cluster,modeling_cluster_custom_dsb_frame)
        uds.initialize_user_defined_dsb()
        self.user_dsb_selector_list.append(uds)


class User_dsb_selector:
    """
    Each modeling cluster will be used to build an object of this class. It will be used to let
    users define custom disulfides bridges in the model chains.
    """
    def __init__(self, pymod_object, modeling_cluster, target_widget):
        self.modeling_cluster = modeling_cluster
        self.target_widget = target_widget
        self.pymod_object = pymod_object


    # Build the initial row in the user-defined disulfide bridges frame.
    def initialize_user_defined_dsb(self):
        # For the rows.
        self.user_disulfides_row_counter = 0
        self.target_list_of_cysteines = []
        # The target chain sequence.
        self.target = self.modeling_cluster.target.my_sequence
        # This is going to contain User_disulfide_combo objects.
        self.list_of_disulfide_combos = []

        # This list is going to contain info about disulfide bridges defined by the user through the
        # GUI. It is going to contain elements like this [[41,xx],[58,yy]] (the numbers are the
        # position of the target cysteines in both the sequence and the alignment).
        self.user_defined_disulfide_bridges = []
        # Builds an interface to let the user define additional dsb only for targets which have at
        # least two CYS residues.
        if self.modeling_cluster.target_with_cys:
            for (k,r) in enumerate(str(self.target).replace("-","")):
                if r == "C":
                    cys = {"position": k + 1,
                        "alignment-position": pmsm.get_residue_id_in_aligned_sequence(self.target,k),
                        "state":"free"}
                    self.target_list_of_cysteines.append(cys)
            self.combobox_frame = Frame(self.target_widget,bg="black")
            self.combobox_frame.grid(row=1)
            # If the target sequence has at least two cys, then creates the comboboxes.
            first = User_disulfide_combo(self.pymod_object,
                    self.user_disulfides_row_counter, self.target_list_of_cysteines,
                    self.combobox_frame, self)
            self.list_of_disulfide_combos.append(first)


    # This is called when the "Add" button to add a user-defined disulfide is pressed.
    def add_new_user_disulfide(self):

        # Checks that both the comboboxes have been used to select a cys.
        if (self.list_of_disulfide_combos[-1].cys1_combobox.get() == "" or self.list_of_disulfide_combos[-1].cys2_combobox.get() == ""):
            txt = "You have to select two cysteines residue to define a disulfide bridge!"
            tkMessageBox.showwarning("Warning", txt,parent=self.pymod_object.modeling_window)

        # Checks that the same cys has not been selected in both comboboxes.
        elif (self.list_of_disulfide_combos[-1].cys1_combobox.get() == self.list_of_disulfide_combos[-1].cys2_combobox.get()):
            txt = "You cannot select the same cysteine to form a disulfide bridge!"
            tkMessageBox.showwarning("Message", txt,parent=self.pymod_object.modeling_window)

        # Checks that the selected cys are not engaged in other bridges.
        # ...

        # If the two cys are free to form a bridge, then adds the new bridge and updates the
        # frame with a new combobox row.
        else:
            self.user_disulfides_row_counter += 1
            # Adds the new row with comboboxes and an "Add" button.
            new_ds_combo = User_disulfide_combo(
                self.pymod_object,
                self.user_disulfides_row_counter,
                self.target_list_of_cysteines,
                self.combobox_frame,
                self)
            # Activates the previous row and returns the name of the 2 selected cys.
            cysteines = self.list_of_disulfide_combos[-1].activate()
            # Finishes and adds the new row.
            self.list_of_disulfide_combos.append(new_ds_combo)
            # Adds the cys pair to the self.user_defined_disulfide_bridges, which is going to be
            # used in the perform_modelization() method.
            self.user_defined_disulfide_bridges.append(cysteines)
            # self.print_user_ds_list()
        self.pymod_object.disulfides_scrolled_frame.reposition()


    # This is called when the "Remove" button is pressed.
    def remove_user_disulfide(self,id_to_remove):
        # Removes the right row of comboboxes.
        for r in self.list_of_disulfide_combos:
            if r.id == id_to_remove:
                # Deactivate and get the right bridge to remove.
                dsb_to_remove = r.deactivate()
                # Finishes to adds the new row.
                self.list_of_disulfide_combos.remove(r)
                # Also removes the bridge from the self.user_defined_disulfide_bridges.
                self.user_defined_disulfide_bridges.remove(dsb_to_remove)
        self.pymod_object.disulfides_scrolled_frame.reposition()


class User_disulfide_combo:
    """
    Class for building in the 'Disulfide' page in the modeling window a "row" with two comboboxes and
    a button to add or remove a user defined disulfide bridge to be included in the model.
    """
    # This is used in the constructor when a new combobox row is created.
    id_counter = 0

    def __init__(self,pymod_object,row,cys_list,target_widget,selector):
        # Selected have the "Add" button, unselected have the "Remove" button.
        self.selected = False
        self.id = User_disulfide_combo.id_counter
        User_disulfide_combo.id_counter += 1
        # Row that is used in the grid method of the widget.
        self.row = row
        # The list of cysteines residues of the target sequence.
        self.cys_list = cys_list
        # The list of strings that is going to appear on the scrollable menus of the comboboxes.
        self.scrollable_cys_list = []
        for cys in self.cys_list:
            self.scrollable_cys_list.append("CYS" + str(cys["position"]))
        self.target_widget = target_widget
        self.pymod_object = pymod_object
        self.selector = selector
        # Creates the first row with two comboboxes.
        self.create_combobox_row()

    def create_combobox_row(self):
        """
        Builds a row with two comboboxes an "Add" button.
        """
        # First CYS combobox.
        self.cys1_combobox = Pmw.ComboBox(
            self.target_widget, label_text = 'Select the first CYS:', labelpos = 'nw',
            selectioncommand = self.select_cys, scrolledlist_items = self.scrollable_cys_list,
            history = 0)
        # Make the combobox entries not editable.
        self.cys1_combobox.component("entryfield").component("entry").configure(
            state='readonly', readonlybackground= "white", fg="black", bg="white")
        self.cys1_combobox.grid(row = self.row,column = 0, padx=(0,0))

        # Second CYS combobox.
        self.cys2_combobox = Pmw.ComboBox(
            self.target_widget, label_text = 'Select the second CYS:', labelpos = 'nw',
            selectioncommand = self.select_cys, scrolledlist_items =self.scrollable_cys_list,
            history = 0)
        self.cys2_combobox.component("entryfield").component("entry").configure(
            state='readonly', readonlybackground= "white", fg="black", bg="white")
        # It must have a little padding to distantiate it from the first combobox and the
        # buttons.
        self.cys2_combobox.grid(row = self.row,column = 1, padx=(15,15))

        self.update_scrollable_cys_list()

        # "Add" button.
        self.new_disulfides_button = Button(self.target_widget, text="Add", command = self.press_add_button, **button_style_2)
        self.new_disulfides_button.grid(row = self.row, column = 2)

        User_disulfide_combo.id_counter += 1

    # This is launched by the combobox when some cys is selected.
    # It should also be used to change the color of the cys according to their state.
    def select_cys(self,i):
        pass

    # Adjust the cysteine list to be displayed on the combobox with the right colors.
    def update_scrollable_cys_list(self):
        # cys = {"position": seq_counter, "alignment-position":k, "state":"free"}
        for (i,cys) in enumerate(self.cys_list):
            if cys["state"] == "free":
                self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
                self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="black")
            elif cys["state"] == "engaged":
                self.cys1_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
                self.cys2_combobox.component("scrolledlist").component("listbox").itemconfig(i,fg="gray")
            # There must also be a condition used to mark cys residues engaged in disulfides
            # present in the templates.

    # Row with two labels and a "Remove button".
    def create_label_row(self):
        # Create dirst CYS label that tells which cys has been selected.
        cys_label_style = {"height" : 1, "background": 'black', "fg":'red', "anchor": "w", "padx": 20, "pady": 7}
        self.cys1_label = Label(self.target_widget, text = self.text1, **cys_label_style)
        self.cys1_label.grid(row = self.row,column = 0)

        # Second CYS label.
        self.cys2_label = Label(self.target_widget, text = self.text2, **cys_label_style)
        self.cys2_label.grid(row = self.row,column = 1)

        # Adds the "Remove" button.
        self.remove_disulfides_button = Button(self.target_widget, text="Remove", command = self.press_remove_button, **button_style_2)
        self.remove_disulfides_button.grid(row = self.row, column = 2,pady=(5,0))


    def press_add_button(self):
        # This is going to call the method below.
        self.selector.add_new_user_disulfide()

    def activate(self):
        """
        This is going to return the information about which cys have been selected when the "Add"
        button is pressed.
        """
        self.selected = True
        # Removes the "Add" button
        self.new_disulfides_button.destroy()
        # Removes both comboboxes, but before get their values.
        self.text1 = self.cys1_combobox.get()
        # print "GET:", self.text1
        self.cys1_combobox.destroy()
        self.text2 = self.cys2_combobox.get()
        self.cys2_combobox.destroy()
        self.create_label_row()
        return self.text1, self.text2

    def press_remove_button(self):
        # This is going to call the method below.
        self.selector.remove_user_disulfide(self.id)

    def deactivate(self):
        """
        This is going to return the information about which bridge has been removed when "Remove"
        button is pressed.
        """
        self.selected = False
        self.cys1_label.destroy()
        self.cys2_label.destroy()
        self.remove_disulfides_button.destroy()
        return self.text1, self.text2


#####################################################################
# Models energy minimization options.                               #
#####################################################################
class Energy_minimization_common:
    restraints_ckb_grid_options = {"sticky":"nw", "padx":10, "pady":2}
    minimization_checkbutton = {"background": widgets_background_color,
                                "foreground": "white", "selectcolor": "red",
                                "highlightbackground": 'black', "font": "comic 8"}

class Energy_minimization_frame(PyMod_frame, Energy_minimization_common):

    def __init__(self, parent = None, borderwidth=1, relief='groove', pady=5, **configs):
        PyMod_frame.__init__(self, parent, borderwidth=borderwidth, relief=relief, pady=pady, **configs)

        #--------------------------
        # Minimization protocols. -
        #--------------------------
        self.title_label = Label(self, text= "Perform Additional Cycles of Energy Minimization", **modeling_window_option_style)
        self.title_label.pack(side="top", anchor="nw")
        self.steepest_descent_frame = Energy_optimizer_frame(self, minimization_algorithm="Steepest Descent", initial_state=1, initial_cycles=80)
        self.steepest_descent_frame.pack(**pack_options_1)
        self.conjugate_gradients_frame = Energy_optimizer_frame(self, minimization_algorithm="Conjugate Gradients", initial_state=1, initial_cycles=20)
        self.conjugate_gradients_frame.pack(**pack_options_1)
        self.quasi_newton_frame = Energy_optimizer_frame(self, minimization_algorithm="Quasi Newton", initial_state=0, initial_cycles=20)
        self.quasi_newton_frame.pack(**pack_options_1)
        self.md_frame = Md_optimizer_frame(self, minimization_algorithm="Molecular Dynamics", initial_state=0, initial_cycles=50)
        self.md_frame.pack(**pack_options_1)

        self.minimization_algorithms_frames = [self.steepest_descent_frame,
                                               self.conjugate_gradients_frame,
                                               self.quasi_newton_frame,
                                               self.md_frame]

        #--------------------------
        # Restranints to satisfy. -
        #--------------------------
        self.parameters_title = Label(self, text= "Select Components to Mimimize", **modeling_window_option_style)
        self.parameters_title.pack(side="top", anchor="nw")
        self.components_frame = PyMod_frame(self)
        self.components_frame.pack(side="top", anchor="nw")

        self.bonds_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Bonds", **self.minimization_checkbutton)
        self.bonds_checkbutton.grid(row=0, column=0, **self.restraints_ckb_grid_options)
        self.angles_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Angles", **self.minimization_checkbutton)
        self.angles_checkbutton.grid(row=1, column=0, **self.restraints_ckb_grid_options)
        self.dihedrals_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Dihedrals", **self.minimization_checkbutton)
        self.dihedrals_checkbutton.grid(row=2, column=0, **self.restraints_ckb_grid_options)
        self.impropers_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Impropers", **self.minimization_checkbutton)
        self.impropers_checkbutton.grid(row=3, column=0, **self.restraints_ckb_grid_options)

        self.lj_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Lennard-Jones", **self.minimization_checkbutton)
        self.lj_checkbutton.grid(row=0, column=1, **self.restraints_ckb_grid_options)
        self.coulomb_checkbutton = Restraint_checkbutton(self.components_frame, value=1, text="Coulomb", **self.minimization_checkbutton)
        self.coulomb_checkbutton.grid(row=1, column=1, **self.restraints_ckb_grid_options)


        self.list_of_parameters_checkbuttons = [self.bonds_checkbutton,
                                                self.angles_checkbutton,
                                                self.dihedrals_checkbutton,
                                                self.impropers_checkbutton,
                                                self.lj_checkbutton,
                                                self.coulomb_checkbutton]

        self.non_bondend_cutoff_rds = Pmw.EntryField(self.components_frame, value = 4.0,
                                                     label_text = "Non Bonded Cutoff (%s):" % (u"\u212B"), labelpos = "w",
                                                     validate = {'validator' : 'real', 'min' : 1.0, 'max' : 250.0})
        self.non_bondend_cutoff_rds.grid(row=2, column=1, sticky="nw", padx=(15,0), pady=2)
        self.non_bondend_cutoff_rds.component("label").configure(**modeling_window_explanation_padless)
        self.non_bondend_cutoff_rds.component("entry").configure(width=4)


class Restraint_checkbutton(Checkbutton):

    def __init__(self, parent, value=1, **configs):
        self.var = IntVar()
        self.var.set(value)
        Checkbutton.__init__(self, parent, variable=self.var, **configs)

    def getvalue(self):
        return self.var.get()


class Energy_optimizer_frame(PyMod_frame, Energy_minimization_common):
    def __init__(self, parent = None, minimization_algorithm=None, additional_information="", initial_state=0, initial_cycles=1, **configs):
        PyMod_frame.__init__(self, parent, **configs)
        self.minimization_algorithm = minimization_algorithm
        self.additional_information = additional_information
        self.use_var = IntVar()
        self.use_var.set(initial_state)
        self.checkbutton = Checkbutton(self, text="Do ", variable=self.use_var, **self.minimization_checkbutton)
        self.checkbutton.pack(side="left")
        self.iterations_enf = Pmw.EntryField(self, value = initial_cycles, labelpos = None,
                                             validate = {'validator' : 'integer', 'min' : 1, 'max' : 1000})
        self.iterations_enf.component("entry").configure(width = 4)
        self.iterations_enf.pack(side="left")
        self.algorithm_label = Label(self, text= " Steps of %s %s" % (self.minimization_algorithm, self.additional_information), **modeling_window_explanation_padless)
        self.algorithm_label.pack(side="left")

    def check_parameters(self):
        if self.iterations_enf.getvalue() != "":
            return True
        else:
            return False


class Md_optimizer_frame(Energy_optimizer_frame):
    def __init__(self, parent = None, minimization_algorithm=None, additional_information="at temparature (K): ", initial_state=0, initial_cycles=1, **configs):
        Energy_optimizer_frame.__init__(self, parent = parent, minimization_algorithm=minimization_algorithm, additional_information=additional_information, initial_state=initial_state, initial_cycles=initial_cycles, **configs)
        self.temperature_enf = Pmw.EntryField(self, value = 293, labelpos = None,
                                              validate = {'validator' : 'integer', 'min' : 1, 'max' : 1000})
        self.temperature_enf.component("entry").configure(width = 4)
        self.temperature_enf.pack(side="left")

    def check_parameters(self):
        if self.iterations_enf.getvalue() != "" and self.temperature_enf.getvalue():
            return True
        else:
            return False


###################################################################################################
# FUNCTIONS USED THROUGHOUT THE MODULE.                                                           #
###################################################################################################

def align_set_of_widgets(widgets_to_align, input_widget_width=10):
        Pmw.alignlabels(widgets_to_align, sticky="nw")
        align_input_widgets_components(widgets_to_align, input_widget_width)


def align_input_widgets_components(widgets_to_align, input_widgets_width):
    """
    Used to force to the same width all the input components of a list of widgets to align.
    It will be generally used along (and after) the label componets are aligned with
    Pmw.alignlabels().
    """
    map(lambda w: w.set_input_widget_width(input_widgets_width), widgets_to_align)


def get_parent_window(target_widget):
    """
    Returns the parent window Tkinter object of a target widget. Useful when specifiying the parents
    windows in Tkinter dialogs.
    """
    parent_window_name = target_widget.winfo_parent()
    parent_window = target_widget.nametowidget(parent_window_name) # also: _nametowidget
    return parent_window
