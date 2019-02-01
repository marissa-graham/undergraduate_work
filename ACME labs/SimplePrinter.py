from PySide import QtGui, QtCore
import sys

class Printer(QtGui.QWidget):
    def __init__(self):
        super(Printer, self).__init__()
        self._initUI()
        
    def _initUI(self):
        self.textBar = QtGui.QLineEdit()
        self.label = QtGui.QLabel()
        
        self.textBar.returnPressed.connect(self.updateText)
        
        vbox = QtGui.QVBoxLayout()
        
        vbox.addWidget(self.textBar, 0, 0)
        vbox.addWidget(self.label, 1, 0)
        
        self.setLayout(vbox)
        
        self.setGeometry(50, 50, 200, 200)
        self.setWindowTitle("Simple Printer")
        self.show()
        
    def updateText(self):
        self.label.setText(self.textBar.displayText())
        self.textBar.clear()

def main():
    app = QtGui.QApplication(sys.argv)
    p = Printer()
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()