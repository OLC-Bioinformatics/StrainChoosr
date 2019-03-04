#!/usr/bin/env python

from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow, QPushButton, QErrorMessage, QLabel, QSpinBox, \
    QColorDialog, QProgressBar, QRadioButton, QListWidget
from PyQt5.QtGui import QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from strainchoosr import strainchoosr
import subprocess
import tempfile
import shutil
import ete3
import sys
import os

# TODO: Have image that gets created be resizeable/zoomable. Pretty useless on bigger trees.
# TODO: Allow user to pick starting strains.
# TODO: Allow user to select tree weights.


class StrainChoosrThread(QThread):
    # All I know about QThreads comes from here: https://kushaldas.in/posts/pyqt5-thread-example.html
    signal = pyqtSignal('PyQt_PyObject')

    def __init__(self, tree_file, num_strains, tmpdir, progress_bar, orientaion, color):
        QThread.__init__(self)
        self.tree_file = tree_file
        self.num_strains = num_strains
        self.tmpdir = tmpdir
        self.progress_bar = progress_bar
        self.orientation = orientaion
        self.color = color

    def run(self):
        tree = ete3.Tree(self.tree_file)
        self.progress_bar.setValue(0)
        diverse_strains = list()
        self.progress_bar.setValue(1)
        diverse_strains = strainchoosr.find_starting_leaves(tree, diverse_strains)

        while len(diverse_strains) < self.num_strains:
            next_leaf = strainchoosr.find_next_leaf(diverse_strains, tree)
            diverse_strains.append(next_leaf)
            self.progress_bar.setValue(100 * len(diverse_strains)/self.num_strains)

        # Future person looking at this - you may be wondering why launching a subprocess here is necessary at all.
        # Here's why - the underlying ete3 code that renders the tree to image uses PyQt and somewhere in there
        # another PyQt application is launched. Then, when this GUI gets closed, the GUI process is still running
        # and has to be manually killed (even Ctrl+C doesn't work), presumably because only one of the two applications gets closed.
        # I couldn't find a way to kill the ete3 PyQt app via the code, so my hacky solution is to run tree rendering via a subprocess
        # so that my GUI doesn't know anything about the ete3 GUI and therefore whatever interaction was occurring
        # can no longer occur.
        cmd = 'strainchoosr_drawimage {} {} {} "{}" {}'.format(self.tree_file, self.num_strains, self.tmpdir, self.color, self.orientation)
        subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        self.progress_bar.setValue(100)
        self.signal.emit(strainchoosr.get_leaf_names_from_nodes(diverse_strains))


class StrainChoosrGUI(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'StrainChoosr GUI'
        self.left =10
        self.top = 10
        self.width = 800
        self.height = 600
        self.fasta_files = list()
        self.tmpdir = tempfile.mkdtemp()
        if not os.path.exists(self.tmpdir):
            os.makedirs(self.tmpdir)
        self.fasta_button = None
        self.newick_button = None
        self.newick_tree_label = None
        self.strainchoosr_button = None
        self.strain_number_input = QSpinBox(self)
        self.newick_tree = None
        self.color_dialog_button = None
        self.color_label = None
        self.color_label_palette = None
        self.color = 'red'
        self.export_chosen_strains = None
        self.chosen_strains = list()
        self.progress = QProgressBar(self)
        self.file_save_button = None
        self.tree_orientation_rect = None
        self.tree_orientation_circ = None
        self.tree_orient = 'r'
        self.st_thread = StrainChoosrThread(tree_file=None,
                                            num_strains=None,
                                            tmpdir=self.tmpdir,
                                            progress_bar=self.progress,
                                            color=self.color,
                                            orientaion=self.tree_orient)
        self.st_thread.signal.connect(self.st_finished)
        self.save_image_button = None
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        # self.fasta_button = QPushButton('Select FASTA files', self)
        # self.fasta_button.move(100, 100)
        # self.fasta_button.resize(200, 40)
        # self.fasta_button.clicked.connect(self.get_fasta_files)
        self.newick_button = QPushButton('Select Newick Tree', self)
        self.newick_button.resize(200, 40)
        self.newick_tree_label = QLabel(self)
        self.newick_tree_label.move(210, 10)
        self.newick_tree_label.resize(100, 20)
        self.newick_tree_label.setText(str(self.newick_tree))
        self.newick_tree_label.show()
        self.newick_button.clicked.connect(self.get_newick_tree)
        self.strainchoosr_button = QPushButton('Run StrainChoosr', self)
        self.strainchoosr_button.setEnabled(False)
        self.strainchoosr_button.clicked.connect(self.run_strainchoosr)
        self.strainchoosr_button.resize(200, 40)
        self.strainchoosr_button.move(0, 200)
        self.color_dialog_button = QPushButton('Choose Strain Color', self)
        self.color_dialog_button.setEnabled(False)
        self.color_dialog_button.move(0, 245)
        self.color_dialog_button.resize(150, 30)
        self.color_dialog_button.clicked.connect(self.choose_color)
        self.strain_number_input.setEnabled(False)
        self.color_label = QLabel(self)
        self.color_label.move(155, 245)
        self.color_label.setAutoFillBackground(True)
        self.color_label_palette = self.color_label.palette()
        self.color_label_palette.setColor(QPalette.Background, QColor(self.color))
        self.color_label.setPalette(self.color_label_palette)
        self.file_save_button = QPushButton('Save Strains', self)
        self.file_save_button.setEnabled(False)
        self.file_save_button.clicked.connect(self.file_save)
        self.file_save_button.resize(200, 40)
        self.file_save_button.move(400, 450)
        self.save_image_button = QPushButton('Save Image', self)
        self.save_image_button.clicked.connect(self.save_image)
        self.save_image_button.resize(200, 40)
        self.save_image_button.move(600, 450)
        self.save_image_button.setEnabled(False)
        self.progress.move(0, 300)
        self.tree_orientation_rect = QRadioButton('Rectangular', self)
        self.tree_orientation_circ = QRadioButton('Circular', self)
        self.tree_orientation_rect.setChecked(True)
        self.tree_orientation_rect.toggled.connect(lambda: self.orientation_btn_state(self.tree_orientation_rect))
        self.tree_orientation_circ.toggled.connect(lambda: self.orientation_btn_state(self.tree_orientation_circ))
        self.tree_orientation_rect.move(0, 350)
        self.tree_orientation_circ.move(0, 375)
        self.tree_orientation_rect.show()
        self.tree_orientation_circ.show()

    def save_image(self):
        options = QFileDialog.Options()
        filename = QFileDialog.getSaveFileName(self,
                                               'Select Save File',
                                               '',
                                               'Image Files (*.png)',
                                               options=options)[0]
        if filename.endswith('.png'):
            shutil.copy(os.path.join(self.tmpdir, 'image.png'), filename)
        else:
            shutil.copy(os.path.join(self.tmpdir, 'image.png'), filename + '.png')

    def orientation_btn_state(self, b):
        if b.text() == 'Rectangular':
            if b.isChecked():
                self.tree_orient = 'r'
        if b.text() == 'Circular':
            if b.isChecked():
                self.tree_orient = 'c'

    def file_save(self):
        output_file_name = QFileDialog.getSaveFileName(self, 'Save File')[0]
        with open(output_file_name, 'w') as f:
            for strain_name in self.chosen_strains:
                f.write('{strain_name}\n'.format(strain_name=strain_name))

    def choose_color(self):
        color = QColorDialog.getColor()
        self.color = color.name()
        self.color_label_palette.setColor(QPalette.Background, QColor(self.color))
        self.color_label.setPalette(self.color_label_palette)

    def get_fasta_files(self):
        options = QFileDialog.Options()
        self.fasta_files = QFileDialog.getOpenFileNames(self,
                                                        'Select FASTA files',
                                                        '',
                                                        'Fasta Files (*.fasta)',
                                                        options=options)[0]

    def get_newick_tree(self):
        options = QFileDialog.Options()
        file_input = QFileDialog.getOpenFileName(self,
                                                 'Select Newick Tree',
                                                 '',
                                                 'Tree Files (*.nwk)',
                                                 options=options)
        self.newick_tree = file_input[0]
        number_leaves = len(ete3.Tree(self.newick_tree).get_leaves())
        strain_number_label = QLabel(self)
        strain_number_label.setText('Select number of strains to choose')
        strain_number_label.move(110, 150)
        strain_number_label.resize(300, 20)
        strain_number_label.show()
        self.strain_number_input.move(0, 150)
        self.strain_number_input.setMaximum(number_leaves)
        self.strain_number_input.setMinimum(2)
        self.strain_number_input.setEnabled(True)
        self.strainchoosr_button.setEnabled(True)
        self.color_dialog_button.setEnabled(True)
        self.color_dialog_button.show()
        self.color_label.show()
        self.newick_tree_label.setText(os.path.split(self.newick_tree)[1])

    def run_strainchoosr(self):
        if self.newick_tree is None:
            msg = QErrorMessage()
            msg.showMessage('You have not selected a tree, do that first!')
            msg.exec_()
        else:
            self.st_thread.tree_file = self.newick_tree
            self.st_thread.num_strains = self.strain_number_input.value()
            self.st_thread.color = self.color
            self.st_thread.orientation = self.tree_orient
            self.st_thread.start()
            self.strainchoosr_button.setEnabled(False)

    def st_finished(self, result):
        self.strainchoosr_button.setEnabled(True)
        pic = QLabel(self)
        pixmap = QPixmap(os.path.join(self.tmpdir, 'image.png'))
        scaled_pixmap = pixmap.scaled(400, 400, Qt.KeepAspectRatio, Qt.FastTransformation)
        pic.setPixmap(scaled_pixmap)
        pic.resize(scaled_pixmap.width(), scaled_pixmap.height())
        pic.move(400, 10)
        pic.show()
        self.file_save_button.setEnabled(True)
        self.save_image_button.setEnabled(True)

    def closeEvent(self, event):
        shutil.rmtree(self.tmpdir)
        self.close()


def draw_image_wrapper():
    tree_file = sys.argv[1]
    num_strains = int(sys.argv[2])
    output_dir = sys.argv[3]
    color = sys.argv[4]
    orientation = sys.argv[5]
    tree = ete3.Tree(tree_file)
    diverse_strains = strainchoosr.pd_greedy(tree=tree, number_tips=num_strains, starting_strains=[])
    strainchoosr.create_colored_tree_tip_image(tree_to_draw=tree,
                                               output_file=os.path.join(output_dir, 'image.png'),
                                               representatives=strainchoosr.get_leaf_names_from_nodes(diverse_strains),
                                               mode=orientation,
                                               color=color)


def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    strainchoosr_gui = StrainChoosrGUI()
    strainchoosr_gui.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
