from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_NewWindow(object):
    def setupUi(self, NewWindow):
        NewWindow.setObjectName("Form")
        NewWindow.resize(479, 391)
        self.obj_lineEdit = QtWidgets.QLineEdit(NewWindow)
        self.obj_lineEdit.setGeometry(QtCore.QRect(160, 120, 161, 31))
        self.obj_lineEdit.setObjectName("obj_lineEdit")
        self.obj_pushButton = QtWidgets.QPushButton(NewWindow)
        self.obj_pushButton.setGeometry(QtCore.QRect(190, 220, 93, 28))
        self.obj_pushButton.setObjectName("obj_pushButton")

        self.retranslate_Ui(NewWindow)
        QtCore.QMetaObject.connectSlotsByName(NewWindow)

    def retranslate_Ui(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Automatic gene probe"))
        self.obj_lineEdit.setPlaceholderText(_translate("Form", "Sequence Number"))
        self.obj_pushButton.setText(_translate("Form", "Pick"))