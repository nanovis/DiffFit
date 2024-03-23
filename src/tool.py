# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from .tablemodel import TableModel

from chimerax.core.tools import ToolInstance

class TutorialTool(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:user/tools/tutorial.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "DiffFit Viewer"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()        

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QLabel, QLineEdit, QVBoxLayout, QTableView        

        layout = QVBoxLayout()
        layout.addWidget(QLabel("Log this text:"))
        self.line_edit = QLineEdit()

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.line_edit.returnPressed.connect(self.return_pressed)
        layout.addWidget(self.line_edit)

        view = QTableView();
        view.resize(800, 500)
        view.horizontalHeader().setStretchLastSection(True)
        view.setAlternatingRowColors(True)
        view.setSelectionBehavior(QTableView.SelectRows)
        view.clicked.connect(self.table_row_clicked)        

        # placeholder for actual Roden's data
        data = []
        rows = 4
        
        col = []
        col.append("Id")
        col.append("Position")
        col.append("Rotation")
        data.append(col)
        
        for i in range(rows):
            col = []
            col.append(i + 1)
            col.append([i, i, i])
            col.append(i)
            data.append(col)
        
        data[1][2] = [1, 0, 0, 0]
        data[2][2] = [0, 1, 0, 0]
        data[3][2] = [0, 0, 1, 0]
        data[4][2] = [0, 0, 0, 1]
        
        self._data = data;
        
        model = TableModel(data)
        view.setModel(model)
        view.show()  
        layout.addWidget(view)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def return_pressed(self):
        # The use has pressed the Return key; log the current text as HTML
        from chimerax.core.commands import run
        # ToolInstance has a 'session' attribute...
        run(self.session, "log html %s" % self.line_edit.text())
        
        # If no models were given, use all atomic structures
        #from chimerax.atomic import AtomicStructure
        #structures = self.session.models.list(type=AtomicStructure)
        #num_atoms = 0

        # Loop through structures and print atoms
        #for s in structures:
            # We get the list of atoms and transformed atomic coordinates
            # as arrays so that we can limit the number of accesses to
            # molecular data, which is slower than accessing arrays directly
            #atoms = s.atoms
            #coords = atoms.scene_coords

            #print(s.scene_coords)
            #run(self.session, "turn y 1 360 model #1")

            # First line for a structure is the number of atoms
            #run(self.session, "log html atoms=%s" % s.num_atoms)

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.line_edit.clear())
        menu.addAction(clear_action)
        
    def table_row_clicked(self, item):
         position = self._data[item.row() + 1][1];
         rotation = self._data[item.row() + 1][2];
         
         print("setting position {0} and rotation {1}".format(position, rotation))
         
         from chimerax.atomic import AtomicStructure
         structures = self.session.models.list(type=AtomicStructure)
                 
         if len(structures) > 0:         
             s = structures[0]
                          
             #print(dir(s))
             #o = s.scene_position.translation()
             #print(o)
             #s.scene_position.translation()
             
             print(position)
             
             from chimerax.core.commands import run
             #view matrix models #1,0.40692,0.60876,0.68104,1.26,0.29794,0.79324,-0.53104,1.11,-0.56351,0.013184,0.50416,2.549
             run(self.session, "view matrix models #1,0.40692,0.60876,0.68104,{0},0.29794,0.79324,-0.53104,{1},-0.56351,0.013184,0.50416,{2}".format(position[0], position[1], position[2]))
            
            # -----------------------------------------------------------------------------
    
