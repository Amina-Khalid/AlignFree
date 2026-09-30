"""
AlignFree.py
============
AlignFree - Alignment-Free Sequence Comparison Tool

Initial version: pairwise and multiple sequence comparison using k-mer
frequency profiles and Euclidean distance, with UPGMA/Neighbor-Joining
phylogenetic tree construction.

Run:
    python AlignFree.py
(requires PyQt5, biopython, matplotlib, scipy - see requirements.txt)
"""

import sys
from itertools import combinations

from scipy.spatial.distance import euclidean

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QGroupBox, QFormLayout, QLineEdit, QTextEdit, QSpinBox, QPushButton,
    QLabel, QTableWidget, QTableWidgetItem, QMessageBox, QHeaderView,
    QScrollArea
)
from PyQt5.QtCore import Qt


# =============================================================================
# CORE LOGIC
# =============================================================================

def get_kmers(sequence, k):
    """Slice a sequence into overlapping k-mers (substrings of length k)."""
    if not sequence:
        raise ValueError("Sequence is empty.")
    if k <= 0:
        raise ValueError("k-mer length must be a positive integer.")

    sequence = sequence.strip().upper().replace(" ", "").replace("\n", "")

    if len(sequence) < k:
        raise ValueError(
            f"Sequence length ({len(sequence)}) is shorter than k-mer length ({k})."
        )

    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def build_vocabulary(list_of_kmer_lists):
    """Union vocabulary (W_k) from several k-mer lists, sorted for a deterministic order."""
    vocab = set()
    for kmers in list_of_kmer_lists:
        vocab.update(kmers)
    return sorted(vocab)


def count_vector(kmers, vocabulary):
    """Count how many times each vocabulary word appears in `kmers`."""
    counts = {word: 0 for word in vocabulary}
    for kmer in kmers:
        if kmer in counts:
            counts[kmer] += 1
    return [counts[word] for word in vocabulary]


def euclidean_distance(vector_x, vector_y):
    """Euclidean distance between two equal-length numeric vectors."""
    return euclidean(vector_x, vector_y)


def pairwise_comparison(seq_x, seq_y, k):
    """Run the full alignment-free pipeline for two sequences."""
    kmers_x = get_kmers(seq_x, k)
    kmers_y = get_kmers(seq_y, k)

    vocabulary = build_vocabulary([kmers_x, kmers_y])

    vector_x = count_vector(kmers_x, vocabulary)
    vector_y = count_vector(kmers_y, vocabulary)

    distance = euclidean_distance(vector_x, vector_y)

    return {
        "kmers_x": kmers_x,
        "kmers_y": kmers_y,
        "vocabulary": vocabulary,
        "vector_x": vector_x,
        "vector_y": vector_y,
        "distance": distance,
    }


def build_distance_matrix(names, sequences, k):
    """Build an all-vs-all Euclidean distance matrix for a set of sequences."""
    if len(names) != len(sequences):
        raise ValueError("Number of names and sequences must match.")
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required for a distance matrix.")

    all_kmers = [get_kmers(seq, k) for seq in sequences]
    vocabulary = build_vocabulary(all_kmers)
    vectors = [count_vector(kmers, vocabulary) for kmers in all_kmers]

    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        d = euclidean_distance(vectors[i], vectors[j])
        matrix[i][j] = d
        matrix[j][i] = d

    return matrix, vocabulary, vectors


def to_lower_triangle(matrix):
    """Convert a full square matrix into Bio.Phylo's lower-triangular list-of-lists format."""
    lower = []
    for i in range(len(matrix)):
        lower.append([matrix[i][j] for j in range(i + 1)])
    return lower


def build_tree(names, matrix, method="upgma"):
    """Build a phylogenetic tree from a distance matrix using Biopython."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    dm = DistanceMatrix(names=list(names), matrix=to_lower_triangle(matrix))
    constructor = DistanceTreeConstructor()

    method = method.lower()
    if method == "upgma":
        tree = constructor.upgma(dm)
    elif method == "nj":
        tree = constructor.nj(dm)
    else:
        raise ValueError("method must be 'upgma' or 'nj'")

    return tree


# =============================================================================
# GUI: Pairwise Comparison tab
# =============================================================================
class PairwiseTab(QWidget):
    def __init__(self):
        super().__init__()
        self._build_ui()

    def _build_ui(self):
        main_layout = QVBoxLayout(self)

        title = QLabel("Pairwise Alignment-Free Sequence Comparison")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        main_layout.addWidget(title)

        inputs_layout = QHBoxLayout()
        self.seq1_name, self.seq1_text = self._make_sequence_box("Sequence 1", inputs_layout, default_name="X")
        self.seq2_name, self.seq2_text = self._make_sequence_box("Sequence 2", inputs_layout, default_name="Y")
        main_layout.addLayout(inputs_layout)

        controls_layout = QHBoxLayout()
        controls_layout.addWidget(QLabel("k-mer length (k):"))
        self.k_spin = QSpinBox()
        self.k_spin.setMinimum(1)
        self.k_spin.setMaximum(50)
        self.k_spin.setValue(3)
        controls_layout.addWidget(self.k_spin)
        controls_layout.addStretch()

        self.run_button = QPushButton("Calculate Distance")
        self.run_button.clicked.connect(self.run_comparison)
        controls_layout.addWidget(self.run_button)

        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_all)
        controls_layout.addWidget(self.clear_button)

        main_layout.addLayout(controls_layout)

        result_box = QGroupBox("Results")
        result_layout = QVBoxLayout()

        self.distance_label = QLabel("Euclidean Distance: -")
        self.distance_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #1a5276;")
        result_layout.addWidget(self.distance_label)

        result_layout.addWidget(QLabel("K-mer count vectors:"))
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        result_layout.addWidget(self.table)

        self.kmers_label = QLabel("")
        self.kmers_label.setWordWrap(True)
        result_layout.addWidget(self.kmers_label)

        result_box.setLayout(result_layout)
        main_layout.addWidget(result_box)

    def _make_sequence_box(self, title, parent_layout, default_name=""):
        box = QGroupBox(title)
        form = QFormLayout()

        name_edit = QLineEdit()
        name_edit.setText(default_name)
        form.addRow("Gene name:", name_edit)

        seq_edit = QTextEdit()
        seq_edit.setPlaceholderText("Paste gene sequence here, e.g. ATGTGTG")
        seq_edit.setFixedHeight(120)
        form.addRow("Sequence:", seq_edit)

        box.setLayout(form)
        parent_layout.addWidget(box)
        return name_edit, seq_edit

    def run_comparison(self):
        name_x = self.seq1_name.text().strip() or "Sequence1"
        name_y = self.seq2_name.text().strip() or "Sequence2"
        seq_x = self.seq1_text.toPlainText().strip()
        seq_y = self.seq2_text.toPlainText().strip()
        k = self.k_spin.value()

        try:
            result = pairwise_comparison(seq_x, seq_y, k)
        except ValueError as e:
            QMessageBox.warning(self, "Invalid input", str(e))
            return

        self.distance_label.setText(f"Euclidean Distance between {name_x} and {name_y}: {result['distance']:.4f}")

        vocabulary = result["vocabulary"]
        self.table.setRowCount(2)
        self.table.setColumnCount(len(vocabulary))
        self.table.setHorizontalHeaderLabels(vocabulary)
        self.table.setVerticalHeaderLabels([name_x, name_y])

        for col, _ in enumerate(vocabulary):
            self.table.setItem(0, col, QTableWidgetItem(str(result["vector_x"][col])))
            self.table.setItem(1, col, QTableWidgetItem(str(result["vector_y"][col])))

        self.kmers_label.setText(
            f"<b>{name_x} k-mers ({k}):</b> {result['kmers_x']}<br>"
            f"<b>{name_y} k-mers ({k}):</b> {result['kmers_y']}"
        )

    def clear_all(self):
        self.seq1_text.clear()
        self.seq2_text.clear()
        self.distance_label.setText("Euclidean Distance: -")
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.kmers_label.setText("")


# =============================================================================
# GUI: Phylogenetic tree pop-up
# =============================================================================
class TreeDialog(QWidget):
    def __init__(self, tree, title="Phylogenetic Tree", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(800, 600)

        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure
        from Bio import Phylo

        layout = QVBoxLayout(self)
        figure = Figure(figsize=(8, 6))
        canvas = FigureCanvas(figure)
        layout.addWidget(canvas)

        axes = figure.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        canvas.draw()


# =============================================================================
# GUI: Multiple Sequence Comparison tab
# =============================================================================
class MultipleTab(QWidget):
    def __init__(self):
        super().__init__()
        self.sequence_rows = []
        self.last_matrix = None
        self.last_names = None
        self._build_ui()

    def _build_ui(self):
        main_layout = QVBoxLayout(self)

        title = QLabel("Multiple Alignment-Free Sequence Comparison")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        main_layout.addWidget(title)

        top_controls = QHBoxLayout()
        top_controls.addWidget(QLabel("Total number of sequences:"))
        self.n_spin = QSpinBox()
        self.n_spin.setMinimum(2)
        self.n_spin.setMaximum(100)
        self.n_spin.setValue(3)
        top_controls.addWidget(self.n_spin)

        self.build_rows_button = QPushButton("Set")
        self.build_rows_button.clicked.connect(self.build_sequence_rows)
        top_controls.addWidget(self.build_rows_button)
        top_controls.addStretch()
        main_layout.addLayout(top_controls)

        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_content = QWidget()
        self.scroll_layout = QVBoxLayout(self.scroll_content)
        self.scroll_area.setWidget(self.scroll_content)
        main_layout.addWidget(self.scroll_area, stretch=2)

        self.build_sequence_rows()

        run_controls = QHBoxLayout()
        run_controls.addWidget(QLabel("k-mer length (k):"))
        self.k_spin = QSpinBox()
        self.k_spin.setMinimum(1)
        self.k_spin.setMaximum(50)
        self.k_spin.setValue(3)
        run_controls.addWidget(self.k_spin)
        run_controls.addStretch()

        self.run_button = QPushButton("Generate Distance Matrix")
        self.run_button.clicked.connect(self.run_comparison)
        run_controls.addWidget(self.run_button)
        main_layout.addLayout(run_controls)

        matrix_box = QGroupBox("Distance Matrix")
        matrix_layout = QVBoxLayout()
        self.matrix_table = QTableWidget()
        self.matrix_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.matrix_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        matrix_layout.addWidget(self.matrix_table)
        matrix_box.setLayout(matrix_layout)
        main_layout.addWidget(matrix_box, stretch=2)

        tree_controls = QHBoxLayout()
        self.show_tree_button = QPushButton("Show Phylogenetic Tree(s)")
        self.show_tree_button.clicked.connect(self.show_trees)
        self.show_tree_button.setEnabled(False)
        tree_controls.addWidget(self.show_tree_button)
        tree_controls.addStretch()
        main_layout.addLayout(tree_controls)

    def build_sequence_rows(self):
        while self.scroll_layout.count():
            item = self.scroll_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self.sequence_rows = []

        n = self.n_spin.value()
        for i in range(n):
            box = QGroupBox(f"Sequence {i + 1}")
            form = QFormLayout()

            name_edit = QLineEdit()
            name_edit.setText(f"Seq{i + 1}")
            form.addRow("Gene name:", name_edit)

            seq_edit = QTextEdit()
            seq_edit.setFixedHeight(70)
            seq_edit.setPlaceholderText("Paste gene sequence here")
            form.addRow("Sequence:", seq_edit)

            box.setLayout(form)
            self.scroll_layout.addWidget(box)
            self.sequence_rows.append((name_edit, seq_edit))

        self.scroll_layout.addStretch()

    def run_comparison(self):
        names = []
        sequences = []
        for name_edit, seq_edit in self.sequence_rows:
            name = name_edit.text().strip() or f"Seq{len(names) + 1}"
            seq = seq_edit.toPlainText().strip()
            if not seq:
                QMessageBox.warning(self, "Missing sequence", f"Sequence for '{name}' is empty.")
                return
            names.append(name)
            sequences.append(seq)

        k = self.k_spin.value()

        try:
            matrix, vocabulary, vectors = build_distance_matrix(names, sequences, k)
        except ValueError as e:
            QMessageBox.warning(self, "Invalid input", str(e))
            return

        self.last_matrix = matrix
        self.last_names = names

        n = len(names)
        self.matrix_table.setRowCount(n)
        self.matrix_table.setColumnCount(n)
        self.matrix_table.setHorizontalHeaderLabels(names)
        self.matrix_table.setVerticalHeaderLabels(names)

        for i in range(n):
            for j in range(n):
                item = QTableWidgetItem(f"{matrix[i][j]:.4f}")
                item.setTextAlignment(Qt.AlignCenter)
                self.matrix_table.setItem(i, j, item)

        self.show_tree_button.setEnabled(True)

    def show_trees(self):
        if self.last_matrix is None:
            return
        try:
            tree = build_tree(self.last_names, self.last_matrix, method="upgma")
            dialog = TreeDialog(tree, title="Phylogenetic Tree - UPGMA", parent=self)
            dialog.show()
            if not hasattr(self, "_open_windows"):
                self._open_windows = []
            self._open_windows.append(dialog)
        except Exception as e:
            QMessageBox.critical(self, "Tree construction error", str(e))


# =============================================================================
# GUI: Main window
# =============================================================================
class AlignFreeMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AlignFree - Alignment-Free Sequence Comparison Tool")
        self.resize(900, 750)

        header = QLabel("AlignFree — Alignment-Free Sequence Comparison Tool")
        header.setAlignment(Qt.AlignCenter)
        header.setStyleSheet("font-size: 18px; font-weight: bold; color: #154360; padding: 8px;")

        tabs = QTabWidget()
        tabs.addTab(PairwiseTab(), "Pairwise Comparison")
        tabs.addTab(MultipleTab(), "Multiple Sequence Comparison")

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(header)
        layout.addWidget(tabs)

        self.setCentralWidget(central)


def main():
    app = QApplication(sys.argv)
    window = AlignFreeMainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
