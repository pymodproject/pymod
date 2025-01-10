# Copyright 2020 by Giacomo Janson. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Module to display a dialog with a progressbar used when running those functions of PyMod where a
time-consuming process is executed. In this way the GUI will not be blocked during the execution
of the process.
"""

import os
import sys
import subprocess
import time

from pymol.Qt import QtCore, QtWidgets

from pymod_lib.pymod_exceptions import PyModInterruptedProtocol


class PyMod_protocol_thread(QtCore.QThread):
    """
    Class for a 'QThread' to launch a function where a time-consuming process is
    executed. It is used when showing a 'Protocol_exec_dialog' dialog.
    """

    # Signals.
    terminate_thread_signal = QtCore.pyqtSignal(int)
    exception_thread_signal = QtCore.pyqtSignal(Exception)

    def set_params(self, function, args, wait_start, wait_end):
        self.function = function
        self.args = args
        self.wait_start = wait_start
        self.wait_end = wait_end

    def run(self):
        if self.wait_start is not None:
            time.sleep(self.wait_start)

        # Attempts to execute the function in this thread.
        try:
            if type(self.args) is dict:
                self.function(**self.args)
            else:
                self.function(*self.args)
            self._wait_end()
            self.terminate_thread_signal.emit(0) # Terminate sucessully.

        # If there was an error, emit the exception, so that it can be handled in the main thread.
        except Exception as e:
            self._wait_end()
            self.exception_thread_signal.emit(e) # Terminate by raising an exception.

    def _wait_end(self):
        if self.wait_end is not None:
            time.sleep(self.wait_end)


class Protocol_exec_dialog(QtWidgets.QDialog):
    """
    A dialog with a progress bar used when running long and time-consuming functions in PyMod.
    """

    is_pymod_window = True

    def __init__(self, app, pymod, function, args,
                 wait_start=0.15,
                 wait_end=0.15,
                 wait_close=0.20,
                 title="Running a new protocol",
                 label_text="Running. Please wait for the protocol to complete.",
                 lock=False,
                 lock_title="Can not Exit",
                 lock_message="Can not safely exit. Please wait for the threa to complete.",
                 progress=True,
                 stdout_silence=False,
                 stdout_filepath=None,
                 backend="qthread"):

        QtWidgets.QDialog.__init__(self, parent=app)

        self.main_window = app
        self.pymod = pymod
        self.function = function
        self.args = args
        self.wait_start = wait_start # Time to wait before the protocol is actually executed.
        self.wait_end = wait_end # Time to wait after the protocol completes.
        self.wait_close = wait_close # Time to wait before the dialog closes.
        self.title = title
        self.label_text = label_text

        self.lock = lock
        self.lock_title = lock_title
        self.lock_message = lock_message
        self.progress = progress
        self.error = None

        if not backend in ("qthread", "python"):
            raise KeyError("Unknown 'backend': %s" % backend)
        self.backend = backend

        self.initUI()
        self.setModal(True) # Set the type of dialog.

        # Revert the terminal output to a file.
        self.stdout_filepath = stdout_filepath
        # Silence all stdout.
        self.stdout_silence = stdout_silence
        self._revert_stdout = self.stdout_silence or self.stdout_filepath is not None

        # Protocol thread.
        if self.backend == "qthread":
            self.p_thread = PyMod_protocol_thread()
            self.p_thread.set_params(self.function, self.args, self.wait_start, self.wait_end)
            self.p_thread.terminate_thread_signal.connect(self.on_terminate_thread_signal)
            self.p_thread.exception_thread_signal.connect(self.on_exception_thread_signal)
            # self.worker = QtCore.QObject()
            # self.worker.moveToThread(self.p_thread)

        elif self.backend == "python":
            # self.p_thread = threading.Thread(target=self.function,
            #                                  args=self.args,
            #                                  daemon=True)
            raise NotImplementedError


    def initUI(self):

        self.setWindowTitle(self.title)

        vertical_layout = QtWidgets.QVBoxLayout()
        self.thread_progressbar = QtWidgets.QProgressBar(self)
        if self.progress:
            self.thread_progressbar.setMinimum(0)
            self.thread_progressbar.setMaximum(0)
        progressbar_label = "Computing..." # "Wait for the protocol to complete."
        self.thread_progressbar.setFormat(progressbar_label)
        self.thread_progressbar.setValue(0)
        vertical_layout.addWidget(self.thread_progressbar)

        self.thread_progress_label = QtWidgets.QLabel(self.label_text, self)
        self.thread_progress_label.setWordWrap(True)
        vertical_layout.addWidget(self.thread_progress_label)

        # Button for canceling the execution of the thread.
        horizontal_layout = QtWidgets.QHBoxLayout()
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_button_click)
        self.cancel_button.setEnabled(not self.lock)
        horizontal_layout.addWidget(self.cancel_button)

        vertical_layout.addLayout(horizontal_layout)

        self.setLayout(vertical_layout)


    def exec_(self):
        """
        Opens the dialog and performs checks on the errors at the end.
        """
        # Redirect stdout.
        if self._revert_stdout:
            original_stdout = sys.stdout
            if self.stdout_silence:
                sys.stdout = open(os.devnull, "w")
            else:
                sys.stdout = open(self.stdout_filepath, "w")

        # Starts the thread before showing the progress dialog.
        self.p_thread.start()
        QtWidgets.QDialog.exec_(self)

        # Revert stdout.
        if self._revert_stdout:
            sys.stdout.close()
            sys.stdout = original_stdout

        # Complete.
        if self.wait_close is not None:
            time.sleep(self.wait_close)
        self.check_error()


    def check_error(self):
        """
        Checks if there was an exception in the child thread.
        """
        if self.error is None:
            return None
        else: # If there was an exception, then raise it in the main thread.
            raise self.error


    def get_cancel_exception(self):
        """
        Getr the exception raised when interrupting a protocol.
        """
        return PyModInterruptedProtocol("Interrupted the protocol")


    # Interactions with the buttons.
    def on_cancel_button_click(self):
        """
        Stops the protocol thread and provides an exception to stop the execution in the main thread.
        """
        self.on_exception_thread_signal(self.get_cancel_exception())


    def closeEvent(self, evnt):

        # Closes the dialog.
        if not self.lock:
            # The user clicked the "close" button on the window.
            if evnt.spontaneous():
                self._on_exception_thread_signal(self.get_cancel_exception())
            super(Protocol_exec_dialog, self).closeEvent(evnt)

        # Can not close the dialog when the user clicks on it.
        else:
            # Inactivates the "close" button on the window.
            if evnt.spontaneous():
                evnt.ignore()


    # Interactions with the thread.
    def on_terminate_thread_signal(self, status):
        """
        The thread has successfully terminated, closes the dialog.
        """
        if status == 0:
            self.close()
        else:
            raise NotImplementedError


    def on_exception_thread_signal(self, e):
        """
        An exception has been raised while executing the protocol. Stores the exception
        raised in the protocol thread, so that it can be raised in the main thread.
        Then closes the dialog, so that the exception can be raised.
        """
        self._on_exception_thread_signal(e)
        self.close()

    def _on_exception_thread_signal(self, e):
        self.error = e
        self._terminate_threads()

    def _terminate_threads(self):
        if self.p_thread.isRunning():
            # self.worker.stop()
            self.p_thread.terminate()
            # self.p_thread.quit()
            # self.p_thread.wait()


    def keyPressEvent(self, event):
        """
        By overriding this method, the dialog will not close when pressing the "esc" key.
        """
        if event.key() == QtCore.Qt.Key_Escape:
            pass
        else:
            QtWidgets.QDialog.keyPressEvent(self, event)
