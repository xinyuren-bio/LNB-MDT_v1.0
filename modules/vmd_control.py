import sys
import os
import time
import pandas as pd
from socket import socket, AF_INET, SOCK_STREAM
import subprocess
from PySide6.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, \
    QSizePolicy, QFrame, QHBoxLayout, QVBoxLayout
from PySide6.QtCore import Qt
from PySide6.QtGui import QDragEnterEvent, QDropEvent

# VMD 命令类
class VMDCommands:
    @staticmethod
    def gotoFrame(frame):
        return f"animate goto {frame}"

    @staticmethod
    def highlightResid(resids):
        resid_str = " ".join(map(str, resids))
        return f"mol delrep 1 0; mol selection resid {resid_str}; mol representation VDW; mol color ColorID 1; mol addrep 0"

# VMD TCP 控制类
class VMDTcp:
    def __init__(self, rctl_path, vmd_path):
        self.rctl = rctl_path
        self.vmd_path = vmd_path
        self.HOST = 'localhost'
        self.PORT = 5050
        self.ADDR = (self.HOST, self.PORT)
        self.tcpClientSocket = None
        self.vmd_process = None

    def attemptConnection(self):
        max_attempts = 5
        for attempt in range(max_attempts):
            try:
                self.tcpClientSocket = socket(AF_INET, SOCK_STREAM)
                self.tcpClientSocket.connect(self.ADDR)
                return 0
            except ConnectionRefusedError:
                time.sleep(1)
                if attempt == max_attempts - 1:
                    return -1

    def start(self):
        if not os.path.exists(self.rctl):
            raise FileNotFoundError(f"remote_ctl.tcl not found at {self.rctl}")
        if not os.path.exists(self.vmd_path):
            raise FileNotFoundError(f"VMD executable not found at {self.vmd_path}")
        self.vmd_process = subprocess.Popen([self.vmd_path, "-e", self.rctl])
        return self.attemptConnection()

    def send_command(self, cmd):
        if self.tcpClientSocket:
            self.tcpClientSocket.send((cmd + "\n").encode())

    def stop(self):
        if self.tcpClientSocket:
            self.send_command("quit")
            self.tcpClientSocket.close()
        if self.vmd_process:
            self.vmd_process.terminate()

# CSV 读取函数
def read_excel_vmd(file_path):
    try:
        comments = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    comments.append(line.strip()[2:])
                else:
                    break

        df = pd.read_csv(file_path, skiprows=len(comments), header=0)
        valid_comments = comments[1] if len(comments) > 1 else ""
        return valid_comments, df
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return None, None

class VMDControlPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.vmd_process = None
        self.vmd_tcp = None
        self.setupUi()
        self.initConnections()
        self.initState()

    def setupUi(self):
        page_vmd = self.parent()
        if not page_vmd:
            print("Error: Parent widget (page_vmd) not found.")
            return

        self.startButton = page_vmd.findChild(QPushButton, "vmd_btn_start")
        self.stopButton = page_vmd.findChild(QPushButton, "vmd_btn_stop")
        self.table = page_vmd.findChild(QTableWidget, "vmd_tablewidget")
        self.infoLabel = page_vmd.findChild(QLabel, "vmd_label")

        scroll_area = page_vmd.findChild(QWidget, "qt_scrollarea")
        if scroll_area:
            scroll_area.setAcceptDrops(False)
            viewport = scroll_area.findChild(QWidget, "qt_scrollarea_viewport")
            if viewport:
                viewport.setAcceptDrops(False)

        if not all([self.startButton, self.stopButton, self.table, self.infoLabel]):
            print("Error: Some UI elements not found. Check object names in Qt Designer.")
            print("Available children in page_vmd:", [child.objectName() for child in page_vmd.findChildren(QWidget)])
            return

        self.table.setSelectionMode(QTableWidget.ExtendedSelection)
        self.table.setVisible(True)
        self.table.setRowCount(0)
        self.table.setColumnCount(0)

        self.stopButton.setEnabled(False)
        self.infoLabel.setText("Click 'Start VMD' to launch VMD, then drag and drop a CSV file")

        self.setAcceptDrops(True)  # 启用拖放

    def initConnections(self):
        self.startButton.clicked.connect(self.pushStartVMD)
        self.stopButton.clicked.connect(self.pushStopVMD)
        self.table.selectionModel().selectionChanged.connect(self.onSelectionChanged)

    def initState(self):
        self.vmd_running = False
        self.vmd_process = None
        self.vmd_tcp = None
        self.df = None
        self.valid_comments = None

    def pushStartVMD(self):
        try:
            vmd_path = "C:/Program Files/VMD/vmd.exe"
            rctl_path = "path/to/remote_ctl.tcl"  # 替换为 remote_ctl.tcl 的实际路径
            self.vmd_tcp = VMDTcp(rctl_path, vmd_path)
            if self.vmd_tcp.start() == 0:
                self.vmd_running = True
                self.startButton.setEnabled(False)
                self.stopButton.setEnabled(True)
                self.infoLabel.setText("VMD is running. Drag and drop a CSV file to load data.")
                print("VMD started successfully.")
            else:
                self.infoLabel.setText("Failed to connect to VMD.")
                print("Failed to connect to VMD.")
        except Exception as e:
            self.infoLabel.setText(f"Failed to start VMD: {str(e)}")
            print(f"Error starting VMD: {e}")

    def pushStopVMD(self):
        if self.vmd_tcp:
            try:
                self.vmd_tcp.stop()
                self.vmd_tcp = None
                self.vmd_running = False
                self.startButton.setEnabled(True)
                self.stopButton.setEnabled(False)
                self.infoLabel.setText("VMD stopped. Click 'Start VMD' to launch again.")
                print("VMD stopped successfully.")
            except Exception as e:
                self.infoLabel.setText(f"Failed to stop VMD: {str(e)}")
                print(f"Error stopping VMD: {e}")
        else:
            self.infoLabel.setText("VMD is not running.")

    def onSelectionChanged(self, selected, deselected):
        selected_rows = [index.row() for index in self.table.selectionModel().selectedRows()]
        print(f"Selected rows: {selected_rows}")
        if self.vmd_running and self.vmd_tcp and self.df is not None:
            for row in selected_rows:
                frame = int(self.df.iloc[row]["Frame"])
                resids = list(map(int, self.df.iloc[row]["Resid"].split()))
                self.vmd_tcp.send_command(VMDCommands.gotoFrame(frame))
                self.vmd_tcp.send_command(VMDCommands.highlightResid(resids))

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        print("dropEvent triggered in VMDControlPanel")
        if not self.vmd_running:
            self.infoLabel.setText("Please start VMD before dropping a file.")
            event.ignore()
            return

        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return

        file_path = urls[0].toLocalFile()
        print(f"File path: {file_path}")
        if file_path.lower().endswith('.csv'):
            self.loadCSV(file_path)
            event.acceptProposedAction()
        else:
            self.infoLabel.setText("Please drop a CSV file.")
            event.ignore()

    def loadCSV(self, file_path):
        print("loadCSV called")
        try:
            if not os.path.exists(file_path):
                self.infoLabel.setText("CSV file not found!")
                return

            self.valid_comments, self.df = read_excel_vmd(file_path)
            if self.df is None:
                self.infoLabel.setText("Error: Failed to read CSV file.")
                return

            # 调整列名以匹配文件中的大小写
            self.df.rename(columns={'Resid': 'resid', 'Resname': 'resname'}, inplace=True)
            # 忽略 resname 和 coordinations 列
            self.df = self.df.drop(columns=['resname', 'coordinations'], errors='ignore')
            frame_cols = [col for col in self.df.columns if col != 'resid']
            self.displayData(frame_cols)
            self.infoLabel.setText(f"CSV loaded successfully. Valid comment: {self.valid_comments}")
        except Exception as e:
            self.infoLabel.setText(f"Error loading CSV: {str(e)}")
            print(f"Error loading CSV: {e}")

    def displayData(self, frame_cols):
        print("displayData called")
        self.table.clear()
        self.table.setRowCount(len(self.df))
        self.table.setColumnCount(len(frame_cols) + 1)
        self.table.setHorizontalHeaderLabels(['resid'] + frame_cols)

        for i, row in self.df.iterrows():
            self.table.setItem(i, 0, QTableWidgetItem(str(row['resid'])))
            for j, frame in enumerate(frame_cols):
                self.table.setItem(i, j + 1, QTableWidgetItem(str(row[frame])))

        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.setVisible(True)
        self.table.update()