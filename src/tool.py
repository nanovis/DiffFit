# vim: set expandtab shiftwidth=4 softtabstop=4:
import math
from datetime import datetime

from PyQt6.QtWidgets import QCheckBox
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

from Qt.QtWidgets import QLabel, QPushButton, QLineEdit, QVBoxLayout, QHBoxLayout, QGridLayout, QComboBox, QFrame
from Qt.QtWidgets import QTableView, QSlider, QTabWidget, QGroupBox, QDoubleSpinBox, QSpinBox 
from Qt.QtWidgets import QFileDialog, QSpacerItem
from Qt.QtCore import QSortFilterProxyModel, Qt

from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run
from chimerax.map.volume import volume_list
from chimerax.atomic import AtomicStructure
from chimerax.geometry import Place
from chimerax.ui import MainToolWindow

from chimerax.core.models import Model
from .cluster_viewer import ClusterSphereModel
from chimerax.core.selection import SELECTION_CHANGED
from chimerax.geometry import translation

from .parse_log import look_at_record, look_at_cluster, look_at_MQS_idx, animate_MQS, animate_MQS_2
from .parse_log import simulate_volume, get_transformation_at_record, zero_cluster_density
from .tablemodel import TableModel
from .DiffAtomComp import diff_atom_comp, cluster_and_sort_sqd_fast, diff_fit, conv_volume, numpy2tensor, \
    linear_norm_tensor

import sys
import numpy as np        
import os
import torch
import psutil
import platform
import ast
        

def create_row(parent_layout, left=0, top=0, right=0, bottom=0, spacing=5):
    row_frame = QFrame()
    parent_layout.addWidget(row_frame)
    row_layout = QHBoxLayout(row_frame)
    row_layout.setContentsMargins(left, top, right, bottom)
    row_layout.setSpacing(spacing)
    return row_layout

class DiffFitSettings:    
    def __init__(self):   
        # viewing
        self.view_output_directory: str = "D:\\GIT\\DiffFit\\dev_data\\output\\dev_comp_domain_fit_3_domains_10s20q_new_q"
        self.view_target_vol_path: str = "D:\\GIT\\DiffFit\\dev_data\\input\\domain_fit_demo_3domains\\density2.mrc"
        self.view_structures_directory: str = "D:\\GIT\\DiffFit\dev_data\input\domain_fit_demo_3domains\subunits_cif"
        
        # computing
        self.input_directory: str = "D:\\GIT\\DiffFit\\dev_data\\input\\domain_fit_demo_3domains"
        self.target_vol_path: str = "D:\\GIT\\DiffFit\\dev_data\\input\\domain_fit_demo_3domains\\density2.mrc"
        self.structures_directory: str = "D:\\GIT\\DiffFit\dev_data\input\domain_fit_demo_3domains\subunits_cif"
        self.structures_sim_map_dir: str = "D:\\GIT\\DiffFit\dev_data\input\domain_fit_demo_3domains\subunits_mrc"
        
        self.output_directory: str = "D:\\GIT\\DiffFit\\dev_data\\output"

        self.target_surface_threshold: float = 0.7
        self.min_cluster_size: float = 100
        self.N_shifts: int = 10
        self.N_quaternions: int = 100
        self.negative_space_value: float = -0.5
        self.learning_rate: float = 0.01
        self.N_iters: int = 201        
        self.out_dir_exist_ok: bool = True
        self.conv_loops: int = 10
        self.conv_kernel_sizes: list = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
        self.conv_weights: list = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        
        self.clustering_shift_tolerance : float = 3.0
        self.clustering_angle_tolerance : float = 6.0


class DiffFitTool(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:user/tools/DiffFit.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "DiffFit"
        self.settings = DiffFitSettings()
        
        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.        
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        #self.tool_window.fill_context_menu = self.fill_context_menu

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()

        self.fit_input_mode = "disk file"

        self.fit_result_ready = False
        self.fit_result = None   

        # Register the selection change callback
        self.session.triggers.add_handler(SELECTION_CHANGED, self.selection_callback)
        self.session.triggers.add_handler('graphics update', self.graphics_update_callback)

        self.spheres = None


    def _build_ui(self):
        
        # the base layout is Vertical
        
        tab_widget = QTabWidget()
        tab_widget.setTabPosition(QTabWidget.West)

        # single fit GUI
        single_fit_group = QGroupBox()
        single_fit_group_layout = QVBoxLayout()
        single_fit_group.setLayout(single_fit_group_layout)
        self.build_single_fit_ui(single_fit_group_layout)
        tab_widget.addTab(single_fit_group, "Single")

        # computation GUI
        compute_group = QGroupBox()
        compute_group_layout = QGridLayout()
        compute_group.setLayout(compute_group_layout)
        self.build_compute_ui(compute_group_layout)
        tab_widget.addTab(compute_group, "Compute")

        # device GUI
        device_group = QGroupBox()
        device_group_layout = QVBoxLayout()
        device_group.setLayout(device_group_layout)
        self.build_device_ui(device_group_layout)
        tab_widget.addTab(device_group, "Device")
        self._device_changed()

        # view GUI
        view_group = QGroupBox()
        view_group_layout = QGridLayout()
        view_group.setLayout(view_group_layout)
        self.build_view_ui(view_group_layout)
        tab_widget.addTab(view_group, "View")

        self.tab_widget = tab_widget
        self.tab_view_group = view_group

        # TODO: place where to update the settings
        self.load_settings()
        
        # Set the layout as the contents of our window        
        layout = QVBoxLayout()                
        layout.addWidget(tab_widget)
        self.tool_window.ui_area.setLayout(layout)
        
        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def load_settings(self):
        print("loading settings...")
        
        self.settings.loading = True
    
        #compute
        self.target_vol_path.setText(self.settings.target_vol_path)
        self.structures_dir.setText(self.settings.structures_directory)
        self.structures_sim_map_dir.setText(self.settings.structures_sim_map_dir)
        self.out_dir.setText(self.settings.output_directory)
        self.target_surface_threshold.setValue(self.settings.target_surface_threshold)
        self.min_cluster_size.setValue(self.settings.min_cluster_size)
        self.n_iters.setValue(self.settings.N_iters)
        self.n_shifts.setValue(self.settings.N_shifts)
        self.n_quaternions.setValue(self.settings.N_quaternions)        
        self.negative_space_value.setValue(self.settings.negative_space_value)
        self.learning_rate.setValue(self.settings.learning_rate)
        self.conv_loops.setValue(self.settings.conv_loops)
        self.conv_kernel_sizes.setText("[{0}]".format(','.join(map(str, self.settings.conv_kernel_sizes))))
        self.conv_weights.setText("[{0}]".format(','.join(map(str, self.settings.conv_weights))))
        
        # view
        self.target_vol.setText(self.settings.view_target_vol_path)     
        self.dataset_folder.setText(self.settings.view_output_directory)        
        self.structures_folder.setText(self.settings.view_structures_directory)        
        
        # clustering
        self.clustering_angle_tolerance.setValue(self.settings.clustering_angle_tolerance)
        self.clustering_shift_tolerance.setValue(self.settings.clustering_shift_tolerance)
        
        self.settings.loading = False
    
    def store_settings(self):    
        #print("settings changed...")
        
        if self.settings.loading :
            return
        
        #compute
        self.settings.target_vol_path = self.target_vol_path.text()
        self.settings.structures_directory = self.structures_dir.text()
        self.settings.structures_sim_map_dir = self.structures_sim_map_dir.text()
        self.settings.output_directory = self.out_dir.text()        
        self.settings.target_surface_threshold = self.target_surface_threshold.value()
        self.settings.min_cluster_size = self.min_cluster_size.value()
        self.settings.N_iters = self.n_iters.value()
        self.settings.N_shifts = self.n_shifts.value()
        self.settings.N_quaternions = self.n_quaternions.value()        
        self.settings.negative_space_value = self.negative_space_value.value()
        self.settings.learning_rate = self.learning_rate.value()
        self.settings.conv_loops = self.conv_loops.value()
        
        kernel_sizes = self.conv_kernel_sizes.text().replace(" ", "").strip('][').split(',')
        self.settings.conv_kernel_sizes = [int(s) for s in kernel_sizes]
        
        weights = self.conv_weights.text().replace(" ", "").strip('][').split(',')
        self.settings.conv_weights = [float(s) for s in kernel_sizes]
        
        #view
        self.settings.view_output_directory = self.dataset_folder.text()
        self.settings.view_structures_directory = self.structures_folder.text()
        self.settings.view_target_vol_path = self.target_vol.text()
        
        # clustering
        self.settings.clustering_angle_tolerance = self.clustering_angle_tolerance.value()
        self.settings.clustering_shift_tolerance = self.clustering_shift_tolerance.value()
        
        #print(self.settings)
        #print(self.settings.view_target_vol_path)

    def build_single_fit_ui(self, layout):
        row = QHBoxLayout()
        layout.addLayout(row)

        mol_label = QLabel("Fit")
        row.addWidget(mol_label)

        from chimerax.map import Volume
        from chimerax.atomic import Structure
        from chimerax.ui.widgets import ModelMenuButton
        self._object_menu = om = ModelMenuButton(self.session, class_filter=Structure)
        mlist = self.session.models.list(type=Structure)
        vlist = self.session.models.list(type=Volume)
        if mlist:
            om.value = mlist[0]
        # om.value_changed.connect(self._object_chosen)
        row.addWidget(om)

        iml = QLabel("in map")
        row.addWidget(iml)

        self._map_menu = mm = ModelMenuButton(self.session, class_filter=Volume)
        if vlist:
            mm.value = vlist[0]
        row.addWidget(mm)

        button_fit = QPushButton()
        button_fit.setText("Fit")
        button_fit.clicked.connect(lambda: self.single_fit_button_clicked())
        row.addWidget(button_fit)

        button_options = QPushButton()
        button_options.setText("Options")
        button_options.clicked.connect(lambda: self._show_or_hide_options())
        row.addWidget(button_options)

        row.addStretch()  # This will add a stretchable space to the right

        # Options panel
        row2 = QVBoxLayout()
        layout.addLayout(row2)
        options = self._create_single_fit_options_gui(None)
        row2.addWidget(options)

        layout.addStretch()

    def _single_fit_set_preset(self, shifts=5, quat=20, gaussian=0):
        self._single_fit_n_shifts.setValue(shifts)
        self._single_fit_n_quaternions.setValue(quat)
        self._single_fit_gaussian_loops.setValue(gaussian)

    def _single_fit_result_save_checkbox_clicked(self):
        self._single_fit_out_dir.setEnabled(self._single_fit_result_save_checkbox.isChecked())
        self._single_fit_out_dir_select.setEnabled(self._single_fit_result_save_checkbox.isChecked())

    def _update_smooth_fields(self, value):
        self.smooth_kernel_sizes.setText(str([5] * value))
        self.smooth_weights.setText(str([1.0] * value))

    def _create_single_fit_options_gui(self, parent):

        from chimerax.ui.widgets import CollapsiblePanel
        self._single_fit_options_panel = p = CollapsiblePanel(parent, title=None)
        f = p.content_area


        # resolution row
        row = create_row(f.layout(), top=20)
        single_fit_res_label = QLabel("Use map simulated from atoms, resolution")
        self._single_fit_res = QDoubleSpinBox()
        self._single_fit_res.setValue(5.0)
        self._single_fit_res.setMinimum(0.0001)
        self._single_fit_res.setMaximum(100.0)
        self._single_fit_res.setSingleStep(0.0001)
        self._single_fit_res.setDecimals(4)
        row.addWidget(single_fit_res_label)
        row.addWidget(self._single_fit_res)


        # Preset row
        row = create_row(f.layout(), top=20)
        preset_fast = QPushButton("Fast")
        preset_balanced = QPushButton("Balanced")
        preset_exhaustive = QPushButton("Exhaustive")
        preset_fast.clicked.connect(lambda: self._single_fit_set_preset())
        preset_balanced.clicked.connect(lambda: self._single_fit_set_preset(10, 50, 3))
        preset_exhaustive.clicked.connect(lambda: self._single_fit_set_preset(10, 100, 10))
        row.addWidget(preset_fast)
        row.addWidget(preset_balanced)
        row.addWidget(preset_exhaustive)
        
        
        # Parameter row
        row = create_row(f.layout())
        n_shifts_label = QLabel("# shifts:")
        self._single_fit_n_shifts = QSpinBox()
        self._single_fit_n_shifts.setValue(5)
        self._single_fit_n_shifts.setMinimum(1)
        self._single_fit_n_shifts.setMaximum(500)
        row.addWidget(n_shifts_label)
        row.addWidget(self._single_fit_n_shifts)

        n_quaternions_label = QLabel("# quaternions:")
        self._single_fit_n_quaternions = QSpinBox()
        self._single_fit_n_quaternions.setValue(20)
        self._single_fit_n_quaternions.setMinimum(1)
        self._single_fit_n_quaternions.setMaximum(500)
        row.addWidget(n_quaternions_label)
        row.addWidget(self._single_fit_n_quaternions)
        row.addStretch()


        # Smooth by row
        row = create_row(f.layout())
        smooth_by_label = QLabel("Smooth by:")
        self._smooth_by = QComboBox()
        smooth_methods = ["PyTorch iterative Gaussian",
                          "ChimeraX incremental Gaussian",
                          "ChimeraX iterative Gaussian"]
        self._smooth_by.addItems(smooth_methods)
        row.addWidget(smooth_by_label)
        row.addWidget(self._smooth_by)

        convs_loops_label = QLabel("Smooth loops:")
        self._single_fit_gaussian_loops = QSpinBox()
        self._single_fit_gaussian_loops.setValue(0)
        self._single_fit_gaussian_loops.setMinimum(0)
        self._single_fit_gaussian_loops.setMaximum(50)
        self._single_fit_gaussian_loops.valueChanged.connect(self._update_smooth_fields)
        row.addWidget(convs_loops_label)
        row.addWidget(self._single_fit_gaussian_loops)
        row.addStretch()


        # Smooth kernel size row
        row = create_row(f.layout())
        smooth_kernel_sizes_label = QLabel()
        smooth_kernel_sizes_label.setText("Kernel sizes [list]:")
        self.smooth_kernel_sizes = QLineEdit()
        self.smooth_kernel_sizes.setText("[]")
        row.addWidget(smooth_kernel_sizes_label)
        row.addWidget(self.smooth_kernel_sizes)
        row.addStretch()

        row = create_row(f.layout())
        smooth_weights_label = QLabel()
        smooth_weights_label.setText("Smooth weights [list]:")
        self.smooth_weights = QLineEdit()
        self.smooth_weights.setText("[]")
        row.addWidget(smooth_weights_label)
        row.addWidget(self.smooth_weights)
        row.addStretch()


        # Save result row
        row = create_row(f.layout(), top=20)
        self._single_fit_result_save_checkbox = QCheckBox()
        self._single_fit_result_save_checkbox.clicked.connect(lambda: self._single_fit_result_save_checkbox_clicked())
        save_res_label = QLabel("Save result to")
        row.addWidget(self._single_fit_result_save_checkbox)
        row.addWidget(save_res_label)

        self._single_fit_out_dir = QLineEdit()
        self._single_fit_out_dir.setDisabled(True)
        self._single_fit_out_dir_select = QPushButton("Select")
        self._single_fit_out_dir.setText("DiffFit_out/single_fit")
        self._single_fit_out_dir_select.setDisabled(True)

        row.addWidget(self._single_fit_out_dir)
        row.addWidget(self._single_fit_out_dir_select)


        return p

    def _show_or_hide_options(self):
        self._single_fit_options_panel.toggle_panel_display()


    def build_compute_ui(self, layout):
        row = 0
               
        target_vol_path_label = QLabel()
        target_vol_path_label.setText("Target Volume:")
        self.target_vol_path = QLineEdit()
        self.target_vol_path.textChanged.connect(lambda: self.store_settings())                
        target_vol_path_select = QPushButton("Select")        
        target_vol_path_select.clicked.connect(lambda: self.select_clicked("Target Volume", self.target_vol_path, False, "MRC Files(*.mrc);;MAP Files(*.map)"))        
        layout.addWidget(target_vol_path_label, row, 0)
        layout.addWidget(self.target_vol_path, row, 1)
        layout.addWidget(target_vol_path_select, row, 2)
        row = row + 1
        
        structures_dir_label = QLabel()
        structures_dir_label.setText("Structures Folder:")
        self.structures_dir = QLineEdit()
        self.structures_dir.textChanged.connect(lambda: self.store_settings())    
        structures_dir_select = QPushButton("Select")        
        structures_dir_select.clicked.connect(lambda: self.select_clicked("Structures folder (containing *.cif)", self.structures_dir))
        layout.addWidget(structures_dir_label, row, 0)
        layout.addWidget(self.structures_dir, row, 1)
        layout.addWidget(structures_dir_select, row, 2)
        row = row + 1
        
        structures_sim_map_dir_label = QLabel()
        structures_sim_map_dir_label.setText("Structures sim-map Folder:")
        self.structures_sim_map_dir = QLineEdit()
        self.structures_sim_map_dir.textChanged.connect(lambda: self.store_settings())   
        structures_sim_map_dir_select = QPushButton("Select")        
        structures_sim_map_dir_select.clicked.connect(lambda: self.select_clicked("Structures sim-map Folder", self.structures_sim_map_dir))
        layout.addWidget(structures_sim_map_dir_label, row, 0)
        layout.addWidget(self.structures_sim_map_dir, row, 1)
        layout.addWidget(structures_sim_map_dir_select, row, 2)
        row = row + 1
        
        out_dir_label = QLabel()
        out_dir_label.setText("Output Folder:")
        self.out_dir = QLineEdit()
        self.out_dir.textChanged.connect(lambda: self.store_settings()) 
        out_dir_select = QPushButton("Select")        
        out_dir_select.clicked.connect(lambda: self.select_clicked("Output Folder", self.out_dir))        
        layout.addWidget(out_dir_label, row, 0)
        layout.addWidget(self.out_dir, row, 1)
        layout.addWidget(out_dir_select, row, 2)
        row = row + 1
        
        target_surface_threshold_label = QLabel()
        target_surface_threshold_label.setText("Target surface threshold:")
        self.target_surface_threshold = QDoubleSpinBox()
        self.target_surface_threshold.setMinimum(0.0)
        self.target_surface_threshold.setMaximum(20.0)
        self.target_surface_threshold.setSingleStep(0.01)
        self.target_surface_threshold.setDecimals(4)
        self.target_surface_threshold.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(target_surface_threshold_label, row, 0)
        layout.addWidget(self.target_surface_threshold, row, 1, 1, 2)
        row = row + 1
        
        min_cluster_size_label = QLabel()
        min_cluster_size_label.setText("Min island size:")
        self.min_cluster_size = QSpinBox()
        self.min_cluster_size.setMinimum(1)
        self.min_cluster_size.setMaximum(500)
        self.min_cluster_size.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(min_cluster_size_label, row, 0)
        layout.addWidget(self.min_cluster_size, row, 1, 1, 2)
        row = row + 1
        
        
        # advanced
        # TODO:
        n_iters_label = QLabel()
        n_iters_label.setText("# iters:")
        self.n_iters = QSpinBox()
        self.n_iters.setMinimum(1)
        self.n_iters.setMaximum(500)
        self.n_iters.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(n_iters_label, row, 0)
        layout.addWidget(self.n_iters, row, 1, 1, 2)
        row = row + 1
        
        n_shifts_label = QLabel()
        n_shifts_label.setText("# shifts:")
        self.n_shifts = QSpinBox()
        self.n_shifts.setMinimum(1)
        self.n_shifts.setMaximum(500)
        self.n_shifts.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(n_shifts_label, row, 0)
        layout.addWidget(self.n_shifts, row, 1, 1, 2)
        row = row + 1
        
        n_quaternions_label = QLabel()
        n_quaternions_label.setText("# quaternions:")
        self.n_quaternions = QSpinBox()
        self.n_quaternions.setMinimum(1)
        self.n_quaternions.setMaximum(500)
        self.n_quaternions.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(n_quaternions_label, row, 0)
        layout.addWidget(self.n_quaternions, row, 1, 1, 2)
        row = row + 1
        
        negative_space_value_label = QLabel()
        negative_space_value_label.setText("Negative space value:")
        self.negative_space_value = QDoubleSpinBox()
        self.negative_space_value.setMinimum(-100)
        self.negative_space_value.setMaximum(100)
        self.negative_space_value.setSingleStep(0.1)
        self.negative_space_value.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(negative_space_value_label, row, 0)
        layout.addWidget(self.negative_space_value, row, 1, 1, 2)
        row = row + 1
        
        learning_rate_label = QLabel()
        learning_rate_label.setText("Learning rate:")
        self.learning_rate = QDoubleSpinBox()
        self.learning_rate.setMinimum(0.00001)
        self.learning_rate.setMaximum(10)
        self.learning_rate.setSingleStep(0.01)
        self.learning_rate.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(learning_rate_label, row, 0)
        layout.addWidget(self.learning_rate, row, 1, 1, 2)
        row = row + 1
        
        convs_loops_label = QLabel()
        convs_loops_label.setText("Conv. loops:")
        self.conv_loops = QSpinBox()
        self.conv_loops.setMinimum(1)
        self.conv_loops.setMaximum(500)
        self.conv_loops.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(convs_loops_label, row, 0)
        layout.addWidget(self.conv_loops, row, 1, 1, 2)
        row = row + 1
        
        conv_kernel_sizes_label = QLabel()
        conv_kernel_sizes_label.setText("Conv. kernel sizes [list]:")
        self.conv_kernel_sizes = QLineEdit()
        self.conv_kernel_sizes.setText("[5, 5, 5, 5, 5, 5, 5, 5, 5, 5]")
        self.conv_kernel_sizes.textChanged.connect(lambda: self.store_settings())        
        layout.addWidget(conv_kernel_sizes_label, row, 0)
        layout.addWidget(self.conv_kernel_sizes, row, 1, 1, 2)
        row = row + 1
        
        conv_weights_label = QLabel()
        conv_weights_label.setText("Conv. weights [list]:")
        self.conv_weights = QLineEdit()
        self.conv_weights.setText("[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]")
        self.conv_weights.textChanged.connect(lambda: self.store_settings())        
        layout.addWidget(conv_weights_label, row, 0)
        layout.addWidget(self.conv_weights, row, 1, 1, 2)
        row = row + 1
        
        #conv_weights_label = QCheckBox()
        #conv_weights_label.setText("Visualize results:")
        #self.conv_weights = QLineEdit()
        #self.conv_weights.setText("[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]")
        #self.conv_weights.textChanged.connect(lambda: self.store_settings())        
        #layout.addWidget(conv_weights_label, row, 0)
        #layout.addWidget(self.conv_weights, row, 1, 1, 2)
        #row = row + 1
        
        
        button = QPushButton()
        button.setText("Run!")
        button.clicked.connect(lambda: self.run_button_clicked())        
        layout.addWidget(button, row, 1, 1, 2)

    def build_device_ui(self, layout):
        row = QHBoxLayout()
        layout.addLayout(row)

        device_label = QLabel("Device:")
        self._device = QComboBox()
        devices = []
        if torch.cuda.is_available():
            devices += [f"cuda:{i}" for i in range(torch.cuda.device_count())]
        devices += ["cpu"]
        self._device.addItems(devices)
        self._device.currentIndexChanged.connect(lambda: self._device_changed())

        row.addWidget(device_label)
        row.addWidget(self._device)
        row.addStretch()

        row = QHBoxLayout()
        layout.addLayout(row)

        self._device_info_label = QLabel()
        self._device_info_label.setWordWrap(True)
        row.addWidget(self._device_info_label)

        row.addStretch()

        layout.addStretch()


    def build_view_ui(self, layout):
        row = 0

        view_input_mode_label = QLabel("Input mode:")
        self._view_input_mode = QComboBox()
        self._view_input_mode.addItems(["disk file", "interactive"])
        self._view_input_mode.currentIndexChanged.connect(lambda: self._view_input_mode_changed())
        layout.addWidget(view_input_mode_label, row, 0)
        layout.addWidget(self._view_input_mode, row, 1, 1, 2)
        row = row + 1

        target_vol_label = QLabel("Target Volume:")
        self.target_vol = QLineEdit()        
        self.target_vol.textChanged.connect(lambda: self.store_settings())
        self.target_vol_select = QPushButton("Select")
        self.target_vol_select.clicked.connect(lambda: self.select_clicked("Target Volume", self.target_vol, False, "MRC Files(*.mrc);;MAP Files(*.map)"))
        layout.addWidget(target_vol_label, row, 0)
        layout.addWidget(self.target_vol, row, 1)
        layout.addWidget(self.target_vol_select, row, 2)
        row = row + 1
        
        structures_folder_label = QLabel("Structures Folder:")
        self.structures_folder = QLineEdit()
        self.structures_folder.textChanged.connect(lambda: self.store_settings()) 
        self.structures_folder_select = QPushButton("Select")
        self.structures_folder_select.clicked.connect(lambda: self.select_clicked("Structures folder (containing *.cif)", self.structures_folder))
        layout.addWidget(structures_folder_label, row, 0)
        layout.addWidget(self.structures_folder, row, 1)
        layout.addWidget(self.structures_folder_select, row, 2)
        row = row + 1
        
        # data folder - where the data is stored
        dataset_folder_label = QLabel("Data Folder:")
        self.dataset_folder = QLineEdit()    
        self.dataset_folder.textChanged.connect(lambda: self.store_settings())                
        self.dataset_folder_select = QPushButton("Select")
        self.dataset_folder_select.clicked.connect(lambda: self.select_clicked("Data Folder", self.dataset_folder))
        layout.addWidget(dataset_folder_label, row, 0)
        layout.addWidget(self.dataset_folder, row, 1)
        layout.addWidget(self.dataset_folder_select, row, 2)
        row = row + 1
        
        clustering_shift_tolerance_label = QLabel()
        clustering_shift_tolerance_label.setText("Clustering - Shift Tolerance:")
        self.clustering_shift_tolerance = QDoubleSpinBox()
        self.clustering_shift_tolerance.setMinimum(0.0)
        self.clustering_shift_tolerance.setMaximum(100)
        self.clustering_shift_tolerance.setSingleStep(1.0)
        self.clustering_shift_tolerance.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(clustering_shift_tolerance_label, row, 0)
        layout.addWidget(self.clustering_shift_tolerance, row, 1, 1, 2)
        row = row + 1
                
        clustering_angle_tolerance_label = QLabel()
        clustering_angle_tolerance_label.setText("Clustering - Angle Tolerance:")
        self.clustering_angle_tolerance = QDoubleSpinBox()
        self.clustering_angle_tolerance.setMinimum(0.0)
        self.clustering_angle_tolerance.setMaximum(180)
        self.clustering_angle_tolerance.setSingleStep(1.0)
        self.clustering_angle_tolerance.valueChanged.connect(lambda: self.store_settings())        
        layout.addWidget(clustering_angle_tolerance_label, row, 0)
        layout.addWidget(self.clustering_angle_tolerance, row, 1, 1, 2)
        row = row + 1
        
        # init button                
        button = QPushButton()
        button.setText("Load")
        button.clicked.connect(lambda: self.init_button_clicked())        
        layout.addWidget(button, row, 1, 1, 2)
        row = row + 1        

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        #layout.addWidget(QLabel("Log this text:"))
        #self.line_edit = QLineEdit()
        #self.line_edit.returnPressed.connect(self.return_pressed)
        #layout.addWidget(self.line_edit)

        # table view of all the results
        view = QTableView()
        view.resize(800, 500)
        view.horizontalHeader().setStretchLastSection(True)
        view.setAlternatingRowColors(True)
        view.setSelectionBehavior(QTableView.SelectRows)
        view.clicked.connect(self.table_row_clicked)        
        layout.addWidget(view)        
        self.view = view
        layout.addWidget(view, row, 0, 1, 3)
        row = row + 1        

        # simple label for simple statistics
        stats = QLabel()
        stats.setText("Stats: ")
        self.stats = stats
        layout.addWidget(stats, row, 0, 1, 3)
        row = row + 1        
        
        # button panel                
        simulate_volume_label = QLabel("Resolution:")        
        self.simulate_volume_resolution = QDoubleSpinBox()
        self.simulate_volume_resolution.setMinimum(0.001)
        self.simulate_volume_resolution.setMaximum(100.0)
        self.simulate_volume_resolution.setSingleStep(0.1)
        self.simulate_volume_resolution.setValue(4.0)        
        simulate_volume = QPushButton()
        simulate_volume.setText("Simulate Volume")
        simulate_volume.clicked.connect(self.simulate_volume_clicked)        
        layout.addWidget(simulate_volume_label, row, 0)
        layout.addWidget(self.simulate_volume_resolution, row, 1)
        layout.addWidget(simulate_volume, row, 2)
        row = row + 1

        zero_density_threshold_label = QLabel("Threshold:")
        self.zero_density_threshold = QDoubleSpinBox()
        self.zero_density_threshold.setMinimum(0.0)
        self.zero_density_threshold.setMaximum(100.0)
        self.zero_density_threshold.setSingleStep(0.0001)
        self.zero_density_threshold.setDecimals(4)
        self.zero_density_threshold.setValue(0.0)
        zero_density_button = QPushButton()
        zero_density_button.setText("Zero density")
        zero_density_button.clicked.connect(self.zero_density_button_clicked)
        layout.addWidget(zero_density_threshold_label, row, 0)
        layout.addWidget(self.zero_density_threshold, row, 1)
        layout.addWidget(zero_density_button, row, 2)
        row = row + 1

        # saving currently selected object
        save_label = QLabel("Save:")
        save_structure_button = QPushButton()
        save_structure_button.setText("Structure")
        save_structure_button.clicked.connect(self.save_structure_button_clicked)
        save_working_vol_button = QPushButton()
        save_working_vol_button.setText("Working volume")
        save_working_vol_button.clicked.connect(self.save_working_vol_button_clicked)
        layout.addWidget(save_label, row, 0)
        layout.addWidget(save_working_vol_button, row, 1)
        layout.addWidget(save_structure_button, row, 2)
        row = row + 1        
        
        # slider for animation
        progress_label1 = QLabel()
        progress_label1.setText("Step: ")

        progress = QSlider(Qt.Horizontal)
        progress.setTickPosition(QSlider.TicksBelow)
        progress.setTickInterval(1) 
        progress.setMinimum(1)
        progress.setMaximum(1)
        progress.valueChanged.connect(self.progress_value_changed)        
        self.progress = progress              
        
        progress_label = QLabel()
        progress_label.setText("1/1")
        self.progress_label = progress_label

        layout.addWidget(progress_label1, row, 0)
        layout.addWidget(progress, row, 1)
        layout.addWidget(progress_label, row, 2)
        row = row + 1        
        
        # Enable ClusterSphere
        enable_spheres = QPushButton()
        enable_spheres.setText("Point cloud visualization")
        enable_spheres.clicked.connect(self.enable_spheres_clicked)                        
        layout.addWidget(enable_spheres, row, 0, 1, 3)
        row = row + 1

        # slider for ClusterSphere offset
        debug_toggle = True
        if debug_toggle:
            CS_offset_label = QLabel()
            CS_offset_label.setText("Offset: ")
            layout.addWidget(CS_offset_label, row, 0)

        # CS_offset = QSlider(Qt.Horizontal)
        # spheres_default_offset = 100
        # CS_offset.setValue(spheres_default_offset)
        # CS_offset.setMinimum(0)
        # CS_offset.setMaximum(300)
        # CS_offset.valueChanged.connect(self.CS_offset_changed)
        # self.CS_offset = CS_offset
        #
        # CS_offset_value_label = QLabel()
        # CS_offset_value_label.setText(str(CS_offset.value()))
        # self.CS_offset_value_label = CS_offset_value_label
        #
        #
        # layout.addWidget(CS_offset, row, 1)
        # layout.addWidget(CS_offset_value_label, row, 2)
        # row = row + 1

        # slider for ClusterSphere scale factor
        # CS_scale_label = QLabel()
        # CS_scale_label.setText("Diff Scale: ")
        #
        # CS_scale = QSlider(Qt.Horizontal)
        # spheres_default_scale = 40
        # CS_scale.setValue(spheres_default_scale)
        # CS_scale.setMinimum(20)
        # CS_scale.setMaximum(100)
        # CS_scale.valueChanged.connect(self.CS_scale_changed)
        # self.CS_scale = CS_scale
        #
        # CS_scale_value_label = QLabel()
        # CS_scale_value_label.setText(str(CS_scale.value()))
        # self.CS_scale_value_label = CS_scale_value_label
        #
        # layout.addWidget(CS_scale_label, row, 0)
        # layout.addWidget(CS_scale, row, 1)
        # layout.addWidget(CS_scale_value_label, row, 2)
        # row = row + 1

        #def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        #from Qt.QtGui import QAction
        #clear_action = QAction("Clear", menu)
        #clear_action.triggered.connect(lambda *args: self.init_folder.clear())
        #menu.addAction(clear_action)


    def _device_changed(self):
        device = self._device.currentText()
        info_text = ""

        if device.startswith("cuda"):
            device_index = int(device.split(":")[1])
            device_name = torch.cuda.get_device_name(device_index)
            total_memory = torch.cuda.get_device_properties(device_index).total_memory / 1024 ** 3  # Convert to GB
            allocated_memory = torch.cuda.memory_allocated(device_index) / 1024 ** 3  # Convert to GB
            cached_memory = torch.cuda.memory_reserved(device_index) / 1024 ** 3  # Convert to GB

            info_text = (
                f"Device: {device_name}\n"
                f"Total Memory: {total_memory:.2f} GB\n"
                f"Allocated Memory: {allocated_memory:.2f} GB\n"
                f"Cached Memory: {cached_memory:.2f} GB"
            )
        elif device.startswith("cpu"):
            cpu_name = platform.processor()
            cpu_count = psutil.cpu_count(logical=True)
            cpu_freq = psutil.cpu_freq()
            total_memory = psutil.virtual_memory().total / 1024 ** 3  # Convert to GB
            available_memory = psutil.virtual_memory().available / 1024 ** 3  # Convert to GB
            used_memory = psutil.virtual_memory().used / 1024 ** 3  # Convert to GB

            info_text = (
                f"CPU Name: {cpu_name}\n"
                f"CPU Count (Logical): {cpu_count}\n"
                f"CPU Frequency: {cpu_freq.current:.2f} MHz\n"
                f"Total Memory: {total_memory:.2f} GB\n"
                f"Available Memory: {available_memory:.2f} GB\n"
                f"Used Memory: {used_memory:.2f} GB"
            )

        self._device_info_label.setText(info_text)


    def _view_input_mode_changed(self):
        if self._view_input_mode.currentText() == "interactive":
            self.fit_input_mode = "interactive"
            self.target_vol.setEnabled(False)
            self.target_vol_select.setEnabled(False)
            self.structures_folder.setEnabled(False)
            self.structures_folder_select.setEnabled(False)
            self.dataset_folder.setEnabled(False)
            self.dataset_folder_select.setEnabled(False)
        elif self._view_input_mode.currentText() == "disk file":
            self.fit_input_mode = "disk file"
            self.target_vol.setEnabled(True)
            self.target_vol_select.setEnabled(True)
            self.structures_folder.setEnabled(True)
            self.structures_folder_select.setEnabled(True)
            self.dataset_folder.setEnabled(True)
            self.dataset_folder_select.setEnabled(True)


    def table_row_clicked(self, item):        
    
        if item.row() != -1:
            self.select_table_item(item.row())

            # spheres interaction
            self.activate_sphere(item.row())
                    
        return


    def get_table_item_size(self, cluster_idx):
        return int(self.e_sqd_clusters_ordered[cluster_idx, 3])

    def get_table_item_transformation(self, cluster_idx):
        mol_idx = int(self.e_sqd_clusters_ordered[cluster_idx, 0])
        record_idx = int(self.e_sqd_clusters_ordered[cluster_idx, 1])

        N_iter = len(self.e_sqd_log[mol_idx, record_idx])
        iter_idx = int(self.e_sqd_clusters_ordered[cluster_idx, 2])

        return get_transformation_at_record(self.e_sqd_log, mol_idx, record_idx, iter_idx)

    def select_table_item(self, index):
        proxyIndex = self.proxyModel.index(index, 0)
        sourceIndex = self.proxyModel.mapToSource(proxyIndex)

        self.cluster_idx = sourceIndex.row()        
        self.mol_idx = int(self.e_sqd_clusters_ordered[self.cluster_idx, 0])
        self.record_idx = int(self.e_sqd_clusters_ordered[self.cluster_idx, 1])

        N_iter = len(self.e_sqd_log[self.mol_idx, self.record_idx])
        iter_idx = int(self.e_sqd_clusters_ordered[self.cluster_idx, 2])

        #self.transformation = get_transformation_at_record(self.e_sqd_log, self.mol_idx, self.record_idx, iter_idx)
        self.transformation = self.get_table_item_transformation(self.cluster_idx)

        if self.fit_input_mode == "interactive":
            self.mol = self.fit_mol_list[self.mol_idx]
            self.mol.display = True
            self.mol.scene_position = self.transformation
        elif self.fit_input_mode == "disk file":
            self.mol = look_at_record(self.mol_folder, self.mol_idx, self.transformation, self.session)

        self.session.logger.info(f"Cluster size: {int(self.e_sqd_clusters_ordered[self.cluster_idx, 3])}")
        self.session.logger.info(f"Highest metric reached at iter : {iter_idx}")

        self.progress.setMinimum(1)
        self.progress.setMaximum(N_iter)
        self.progress.setValue(iter_idx + 1)        
        
    def select_clicked(self, text, target, save = False, pattern = "dir"):
        fileName = ""
        ext = ""
        
        if save:
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            fileName, ext = QFileDialog.getSaveFileName(target, text, "", pattern, options = options)
            ext = ext[-4:]
            ext = ext[:3]                
        else:
            if pattern == "dir":
                fileName = QFileDialog.getExistingDirectory(target, text)
            elif len(pattern) > 0 :
                options = QFileDialog.Options()
                options |= QFileDialog.DontUseNativeDialog
                fileName, ext = QFileDialog.getOpenFileName(target, text, "", pattern, options = options)   
                ext = ext[-4:]
                ext = ext[:3]                
                
            if len(fileName) > 0:
                print("settings to GUI")
                print(fileName)
                target.setText(fileName)
            
        return fileName, ext
    
    def show_results(self, e_sqd_log, mol_centers):
        if e_sqd_log is None:
            return

        if self.fit_input_mode == "disk file":
            print("clear the volumes from the scene")

            vlist = volume_list(self.session)
            for v in vlist:
                v.delete()

            print("opening the volume...")
            # print(self.settings)
            print(self.settings.view_target_vol_path)
            self.vol = run(self.session, "open {0}".format(self.settings.view_target_vol_path))[0]

            # TODO: define mol_centers

        elif self.fit_input_mode == "interactive":
            self.vol = self.fit_vol
            self.vol.display = True

        N_mol, N_quat, N_shift, N_iter, N_metric = e_sqd_log.shape
        self.e_sqd_log = e_sqd_log.reshape([N_mol, N_quat * N_shift, N_iter, N_metric])
        self.e_sqd_clusters_ordered = cluster_and_sort_sqd_fast(self.e_sqd_log, mol_centers,
                                                                self.settings.clustering_shift_tolerance,
                                                                self.settings.clustering_angle_tolerance)
        
        self.model = TableModel(self.e_sqd_clusters_ordered, self.e_sqd_log)
        self.proxyModel = QSortFilterProxyModel()
        self.proxyModel.setSourceModel(self.model)
        
        self.view.setModel(self.proxyModel)
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(0, Qt.AscendingOrder)
        self.view.reset()
        self.view.show()  
        
        self.stats.setText("stats: {0} entries".format(self.model.rowCount())) 
        
        print("showing the first cluster")
        self.mol_folder = self.settings.view_structures_directory
        self.cluster_idx = 0
        # look_at_cluster(self.e_sqd_clusters_ordered, self.mol_folder, self.cluster_idx, self.session)

    def _create_volume_conv_list(self, vol, smooth_by, smooth_loops, session, negative_space_value=-0.5):
        # From here on, there are three strategies for utilizing gaussian smooth
        # 1. with increasing sDev on the same input volume
        # 2. with the same sDev iteratively
        # Combine 1 & 2
        # Need to do experiment to see which one is better

        volume_conv_list = [None] * (smooth_loops + 1)
        volume_conv_list[0] = vol.full_matrix()

        if smooth_by == "PyTorch iterative Gaussian":
            volume_conv_list[0], _ = numpy2tensor(volume_conv_list[0], self._device.currentText())
            volume_conv_list[0] = linear_norm_tensor(volume_conv_list[0])
            volume_conv_list = conv_volume(volume_conv_list[0],
                                           self._device.currentText(),
                                           smooth_loops,
                                           ast.literal_eval(self.smooth_kernel_sizes.text()),
                                           negative_space_value=negative_space_value,
                                           kernel_type="Gaussian")
            volume_conv_list = [v.squeeze().detach().cpu().numpy() for v in volume_conv_list]
        elif smooth_by == "ChimeraX incremental Gaussian":
            for conv_idx in range(1, smooth_loops + 1):
                vol_gaussian = run(session, f"volume gaussian #{vol.id[0]} sDev {conv_idx}")

                vol_device, _ = numpy2tensor(vol_gaussian.full_matrix(), self._device.currentText())
                vol_device = linear_norm_tensor(vol_device)

                eligible_volume_tensor = vol_device > 0.0
                vol_device[~eligible_volume_tensor] = negative_space_value

                volume_conv_list[conv_idx] = vol_device.squeeze().detach().cpu().numpy()

                vol_gaussian.delete()
        elif smooth_by == "ChimeraX iterative Gaussian":
            kernel_sizes = ast.literal_eval(self.smooth_kernel_sizes.text())
            vol_current = vol
            for conv_idx in range(1, smooth_loops + 1):
                vol_gaussian = run(session, f"volume gaussian #{vol_current.id[0]} sDev {kernel_sizes[conv_idx - 1]}")

                if conv_idx > 1:
                    vol_current.delete()

                vol_device, _ = numpy2tensor(vol_gaussian.full_matrix(), self._device.currentText())
                vol_device = linear_norm_tensor(vol_device)

                eligible_volume_tensor = vol_device > 0.0
                vol_device[~eligible_volume_tensor] = negative_space_value

                volume_conv_list[conv_idx] = vol_device.squeeze().detach().cpu().numpy()

                vol_current = vol_gaussian

            vol_current.delete()


        return volume_conv_list


    def single_fit_button_clicked(self):
        single_fit_timer_start = datetime.now()

        # Prepare mol anv vol
        mol = self._object_menu.value
        self.fit_mol_list = [mol]

        self.fit_vol = self._map_menu.value
        vol_matrix = self.fit_vol.full_matrix()

        # Copy vol and make it clean after thresholding
        vol_copy = self.fit_vol.writable_copy()
        vol_copy_matrix = vol_copy.data.matrix()
        vol_copy_matrix[vol_copy_matrix < self.fit_vol.maximum_surface_level] = 0
        vol_copy.data.values_changed()

        # Smooth the volume
        smooth_loops = self._single_fit_gaussian_loops.value()
        smooth_by = self._smooth_by.currentText()
        volume_conv_list = self._create_volume_conv_list(vol_copy, smooth_by, smooth_loops, self.session)
        vol_copy.delete()

        # Simulate a map for the mol
        from chimerax.map.molmap import molecule_map
        mol_vol = molecule_map(self.session, mol.atoms, self._single_fit_res.value(), grid_spacing=self.fit_vol.data.step[0])

        # Fit
        timer_start = datetime.now()
        mol_centers, self.fit_result = diff_fit(volume_conv_list,
                                   self.fit_vol.data.step,
                                   self.fit_vol.data.origin,
                                   10,
                                   [mol.atoms.coords],
                                   [(mol_vol.full_matrix(), mol_vol.data.step, mol_vol.data.origin)],
                                   N_shifts=self._single_fit_n_shifts.value(),
                                   N_quaternions=self._single_fit_n_quaternions.value(),
                                   save_results=self._single_fit_result_save_checkbox.isChecked(),
                                   out_dir=self._single_fit_out_dir.text(),
                                   out_dir_exist_ok=True,
                                   device=self._device.currentText()
                                   )
        timer_stop = datetime.now()
        print(f"Fit time elapsed: {timer_stop - timer_start}\n\n")

        mol_vol.delete()

        self._view_input_mode.setCurrentText("interactive")
        self._view_input_mode_changed()
        self.fit_result_ready = True
        self.show_results(self.fit_result, mol_centers)

        self.tab_widget.setCurrentWidget(self.tab_view_group)

        timer_stop = datetime.now()
        print(f"Single fit time elapsed: {timer_stop - single_fit_timer_start}\n\n")

        self.select_table_item(0)


    def run_button_clicked(self):
        #import sys
        #sys.path.append('D:\\GIT\\DiffFit\\src')
        print("Running the computation...")

        #target_vol_path = "D:\\GIT\\DiffFit\dev_data\input\domain_fit_demo_3domains\density2.mrc"
        #output_folder = "D:\\GIT\\DiffFit\dev_data\output"
        
        mol_centers, e_sqd_log = diff_atom_comp(
            target_vol_path=self.settings.target_vol_path,
            target_surface_threshold=self.settings.target_surface_threshold,
            min_cluster_size=self.settings.min_cluster_size,
            structures_dir=self.settings.structures_directory,
            structures_sim_map_dir=self.settings.structures_sim_map_dir,
            N_shifts=self.settings.N_shifts,
            N_quaternions=self.settings.N_quaternions,
            negative_space_value=self.settings.negative_space_value,
            learning_rate=self.settings.learning_rate,
            n_iters=self.settings.N_iters,
            out_dir=self.settings.output_directory,
            out_dir_exist_ok=self.settings.out_dir_exist_ok,
            conv_loops=self.settings.conv_loops,
            conv_kernel_sizes=self.settings.conv_kernel_sizes,
            conv_weights=self.settings.conv_weights,
            device=self._device.currentText()
        )

        # copy the directories
        self.target_vol.setText(self.settings.target_vol_path)     
        self.structures_folder.setText(self.settings.structures_directory)
        self.dataset_folder.setText("{0}".format(self.settings.output_directory))
        #print(self.settings)
        
        # output is tensor
        self.show_results(e_sqd_log.detach().cpu().numpy(), mol_centers)
        self.tab_widget.setCurrentWidget(self.tab_view_group)
        self.select_table_item(0)

    def init_button_clicked(self):
        if self.fit_input_mode == "interactive":
            if self.fit_result_ready:
                self.show_results(self.fit_result)
                return
            else:
                from chimerax.log.cmd import log
                log(self.session, text="Fitting is not performed yet.", error_dialog=True)
                return

        if self.settings is None:
            return
            
        datasetoutput = self.settings.view_output_directory        
        
        if len(datasetoutput) == 0:
            print("Specify the datasetoutput folder first!")
            return
                
        print("loading data...")
        e_sqd_log = np.load("{0}\\e_sqd_log.npy".format(datasetoutput))
        
        #print(e_sqd_log)
        self.show_results(e_sqd_log)

    def save_working_vol_button_clicked(self):
        if not self.vol:
            return

        fileName, ext = self.select_clicked(f"Save {self.vol.name} as", self.view, True, "MRC Density map(*.mrc);;CCP4 density map (*.map)")

        self.save_working_volume(fileName, ext)
        
    def save_structure_button_clicked(self):
        if not self.mol:
            return
        
        fileName, ext = self.select_clicked(f"Save {self.mol.name} as", self.view, True, "CIF Files(*.cif);;PDB Files (*.pdb)")
        
        self.save_structure(fileName, ext)
    
    def save_structure(self, targetpath, ext):
        
        if len(targetpath) > 0 and self.mol:
            run(self.session, "save '{0}.{1}' models #{2}".format(targetpath, ext, self.mol.id[0]))

    def save_working_volume(self, targetpath, ext):

        if len(targetpath) > 0 and self.vol:
            run(self.session, "save '{0}.{1}' models #{2}".format(targetpath, ext, self.vol.id[0]))

    def simulate_volume_clicked(self):
        res = self.simulate_volume_resolution.value()
        if self.fit_input_mode == "disk file":
            self.mol_vol = simulate_volume(self.session, self.vol, self.mol_folder, self.mol_idx, self.transformation,
                                           res)
        elif self.fit_input_mode == "interactive":
            from chimerax.map.molmap import molecule_map
            self.mol_vol = molecule_map(self.session, self.mol.atoms, res, grid_spacing=self.vol.data_origin_and_step()[1][0] / 3)

        return
        
    def zero_density_button_clicked(self):      
        if self.vol is None:
            print("You have to select a row first!")
            return
            
        if self.mol_vol is None:
            print("You have to simulate a volume first!")
            return

        work_vol = zero_cluster_density(self.session, self.mol_vol, self.mol, self.vol, self.cluster_idx,
                                        zero_iter=0, threshold=self.zero_density_threshold.value())
        self.vol = work_vol
        return
    
    def progress_value_changed(self):
        progress = self.progress.value()
        self.progress_label.setText("{0}/{1}".format(progress, self.progress.maximum()))
        
        if self.e_sqd_log is not None and self.mol is not None:
            self.transformation = get_transformation_at_record(self.e_sqd_log, self.mol_idx, self.record_idx, progress - 1)
            self.mol.scene_position = self.transformation
    

    # point cloud visualization
    def CS_offset_changed(self):
        self.CS_offset_value_label.setText(str(self.CS_offset.value()))
        self.update_spheres_pos()

    def CS_scale_changed(self):
        self.CS_scale_value_label.setText(str(self.CS_scale.value()))
        sphere_size = 0.3
        if self.spheres:
            child_models = self.spheres.child_models()
            for i in range(len(child_models)):
                child_models[i].radius = sphere_size * math.pow(child_models[i].hit_number, 20.0/(120 - self.CS_scale.value()))

    def get_model_by_name(self, model_name):
        models = self.session.models.list()
        for model in models:            
            if model.name == model_name:
                return model
                
        return None   

    def focus_table_row(self, idx):
        selection_model = self.view.selectionModel()

        index = self.proxyModel.index(idx, 0)
        selection_model.clear()
        selection_model.select(
            index,
            selection_model.Select
        )

        self.view.scrollTo(index, QTableView.PositionAtCenter)        
        self.view.setFocus()

    def selection_callback(self, trigger, changes):        
        selected_models = self.session.selection.models()
                
        # Print the selected models
        if selected_models:
            for model in selected_models:
                if type(model) is ClusterSphereModel:
                    self.focus_table_row(model.id[1] - 1)
                    self.select_table_item(model.id[1] - 1)

    def update_spheres_pos(self):
        if self.spheres:
            child_models = self.spheres.child_models()
            for i in range(len(child_models)):
                child_models[i].position = translation(child_models[i].original_position +
                                                       self.session.main_view.camera.position.axes()[0] *
                                                       self.CS_offset.value())

    def graphics_update_callback(self, trigger, changes):
        # print(self.session.main_view.camera.view_direction())
        self.update_spheres_pos()
        
    def activate_sphere(self, cluster_idx):

        if self.spheres:
            command = 'select #{0}.{1}'.format(self.spheres.id[0], cluster_idx + 1)
            run(self.session, command)             

    # coloring of the sphere
    def get_sphere_color(self, idx, count):
        # TODO: based on some feature
        value = idx / (count - 1)

        g = 255 * (1 - value)
        r = 255 * value
        b = 0
        a = 255

        return [r, g, b, a]

    def enable_spheres_clicked(self):

        if not self.spheres:
            spheres = Model("clusterSpheres", self.session)
            self.session.models.add([spheres])

            entries_count = self.proxyModel.rowCount()

            # map_x_length = abs(self.vol.xyz_bounds()[1][0] - self.vol.xyz_bounds()[0][0])

            sphere_size = 0.3

            mol_center = self.mol.atoms.coords.mean(axis=0)

            spheres_default_offset = 100
            spheres_default_scale = 40
            for entry_id in range(1, entries_count + 1):
                place = self.get_table_item_transformation(entry_id - 1)
                hit_number = self.get_table_item_size(entry_id - 1)
                original_position = place * mol_center
                place_position = (original_position +
                                  self.session.main_view.camera.position.axes()[0] *
                                  spheres_default_offset)
                color = self.get_sphere_color(entry_id - 1, entries_count)

                spheres.add([ClusterSphereModel(str(entry_id), self.session, color, place_position,
                                                sphere_size * math.pow(hit_number, 20.0/(120 - spheres_default_scale)),
                                                original_position, hit_number)])

            self.spheres = spheres

        return