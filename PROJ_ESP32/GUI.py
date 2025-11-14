import tkinter as tk
from tkinter import ttk
import threading
import struct
import serial
import numpy as np
import queue
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import csv
import os
import datetime



# Global variables
update_time = 100  # Time of update in ms
decimal_round = 5
running = False
flag = 0
buffer_size = 60    # In seconds
# serialPort = serial.Serial("COM4", baudrate=1382400, timeout=2)

# Buffers for graphs
heart_phase_values = np.zeros(buffer_size)
breath_phase_values = np.zeros(buffer_size)
time_values = np.arange(-buffer_size + 1, 1, 1)

def initialize_serial():
    """Initialize the serial port."""
    global serialPort
    try:
        serialPort = serial.Serial(
            port="COM6", baudrate=1382400, bytesize=8, timeout=2, stopbits=serial.STOPBITS_TWO
        )
        print("Serial port initialized successfully.")
    except Exception as e:
        print(f"Error initializing serial port: {e}")

class SerialReaderApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Radar Data Reader")
        self.queue = queue.Queue()

        # Innitialize CSV file
        self.csv_file = None
        self.csv_writer = None
        self.csv_filename = None
        self.previous_timestamp = None


        # Create the main frame with responsive layout
        main_frame = ttk.Frame(root)
        main_frame.grid(row=0, column=0, sticky="nsew")
        main_frame.columnconfigure(0, weight=3)  # Graphs take 3/4 of the width
        main_frame.columnconfigure(1, weight=1)  # Labels take 1/4 of the width
        main_frame.rowconfigure(0, weight=1)

        # Left side: Graphs
        graph_frame = ttk.Frame(main_frame)
        graph_frame.grid(row=0, column=0, sticky="nsew")

        self.buffer_size = int(buffer_size * 1000 / update_time)
        global time_values, heart_phase_values, breath_phase_values
        time_values = np.linspace(-buffer_size, 0, self.buffer_size)
        heart_phase_values = np.zeros(self.buffer_size)
        breath_phase_values = np.zeros(self.buffer_size)

        # Create the graphs
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True)
        self.ax1.set_title("Heart Phase", fontsize=10)
        self.ax1.set_ylim(-8, 8)
        self.ax1.set_xlim(-buffer_size, 0)
        self.ax1.grid(True, linestyle='--', alpha=0.5)
        self.ax2.set_title("Breath Phase", fontsize=10)
        self.ax2.set_ylim(-8, 8)
        self.ax2.set_xlabel("Time (seconds)", fontsize=10)
        self.ax2.set_xlim(-buffer_size, 0)
        self.ax2.grid(True, linestyle='--', alpha=0.5)
        self.line1, = self.ax1.plot(time_values, heart_phase_values, "r-", linewidth=1.5)
        self.line2, = self.ax2.plot(time_values, breath_phase_values, "b-", linewidth=1.5)
        self.canvas = FigureCanvasTkAgg(self.fig, master=graph_frame)
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)

        # Right side: Labels (Squares)
        square_frame = ttk.Frame(main_frame)
        square_frame.grid(row=0, column=1, sticky="nsew", rowspan=6, padx=10)
        square_frame.columnconfigure(0, weight=1)

        self.square_labels = []
        for i, label_text in enumerate(["Heart:", "Breath:", "Total:", "Distance:"]):
            label = tk.Label(
                square_frame,
                bg="#f0f0f0",
                font=("Arial", 14, "bold"),
                fg="black",
                borderwidth=2,
                relief="solid",
                width=20,
                height=3,
                anchor="center"
            )
            label.config(text=f"{label_text}\n0.00000")
            label.grid(row=i, column=0, pady=5, sticky="ew")
            self.square_labels.append(label)

        # Add buttons under the labels
        square_frame.rowconfigure(4, weight=1)
        square_frame.rowconfigure(5, weight=1)

        # Start Button (Row 4)
        self.start_button = tk.Button(
            square_frame,
            text="Start Serial",
            command=self.start_plotting,
            bg="yellowgreen",
            activebackground="lightgreen",
            font=("Arial", 14, "bold"),
            fg="black",
            borderwidth=2,
            relief="solid",
            width=20,
            height=3
        )
        
        self.start_button.grid(row=4, column=0, pady=5, sticky="ew")

        # Stop Button (Row 5)
        self.stop_button = tk.Button(
            square_frame,
            text="Stop Serial",
            command=self.stop_plotting,
            bg="tomato",
            activebackground="coral",
            font=("Arial", 14, "bold"),
            fg="black",
            borderwidth=2,
            relief="solid",
            width=20,
            height=3
        )
        self.stop_button.grid(row=5, column=0, pady=5, sticky="ew")
        self.start_button.config(state=tk.NORMAL)
        self.stop_button.config(state=tk.DISABLED)

        # Log area with scrollbar and auto-scroll
        log_frame = ttk.Frame(root)
        log_frame.grid(row=1, column=0, sticky="nsew")
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.text_log = tk.Text(log_frame, height=6, state=tk.DISABLED, wrap=tk.NONE)
        self.text_log.grid(row=0, column=0, sticky="nsew")
        scrollbar_x = ttk.Scrollbar(log_frame, orient=tk.HORIZONTAL, command=self.text_log.xview)
        scrollbar_x.grid(row=1, column=0, sticky="ew")
        scrollbar_y = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.text_log.yview)
        scrollbar_y.grid(row=0, column=1, sticky="ns")
        self.text_log.config(xscrollcommand=scrollbar_x.set, yscrollcommand=scrollbar_y.set)

        # Start the serial reading thread
        self.thread = threading.Thread(target=self.read_serial_thread, daemon=True)
        self.thread.start()

        # Start the UI update loop
        self.update_counter = 0  # Counter to control graph redraw frequency
        self.update_ui()

    def start_plotting(self):
        """Start the serial reading process."""
        global running
        if not running:
            running = True
            self.start_button.config(state=tk.DISABLED)
            self.stop_button.config(state=tk.NORMAL)
            print("Starting serial reading...")

            # Create a new CSV file
            now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            self.csv_filename = f"radar_data_{now}.csv"
            self.csv_file = open(self.csv_filename, mode="w", newline="")
            self.csv_writer = csv.writer(self.csv_file)
            self.csv_writer.writerow(["Timestamp", "HeartPhase", "BreathPhase", "TotalPhase", "Distance"])
            self.current_data = {}

            self.update_ui()

            self.update_ui()

    def stop_plotting(self):
        global running
        if running:
            running = False
            self.stop_button.config(state=tk.DISABLED)
            self.start_button.config(state=tk.NORMAL)
            print("Stopping serial reading...")

             # Close the CSV file
            if self.csv_file:
                self.csv_file.close()
                print(f"Data saved to {self.csv_filename}")
            
    def read_serial_thread(self):
        """Read data from the serial port in a separate thread."""
        global flag, heart_phase_values, breath_phase_values, running
        while True:
            try:
                byte = serialPort.read(1)
                if byte == b'\x0A':
                    next_byte = serialPort.read(1)
                    if next_byte == b'\x14':
                        flag = 0
                    elif next_byte == b'\x13':
                        data = serialPort.read(13)
                        TotalPhase = struct.unpack('<f', data[1:5])[0]
                        BreathPhase = struct.unpack('<f', data[5:9])[0]
                        HeartPhase = struct.unpack('<f', data[9:13])[0]
                        if running:                            
                            self.current_data["timestamp"] = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]
                            self.current_data["HeartPhase"] = HeartPhase
                            self.current_data["BreathPhase"] = BreathPhase
                            self.current_data["TotalPhase"] = TotalPhase
                            self.queue.put({
                                "type": "data",
                                "message": f"Breath: {BreathPhase}, Heart: {HeartPhase}, Total: {TotalPhase}\n",
                                "HeartPhase": HeartPhase,
                                "BreathPhase": BreathPhase,
                                "TotalPhase": TotalPhase
                            })

                            # Update the buffers for the graphs
                            heart_phase_values[:-1] = heart_phase_values[1:]
                            heart_phase_values[-1] = HeartPhase
                            breath_phase_values[:-1] = breath_phase_values[1:]
                            breath_phase_values[-1] = BreathPhase
                        flag = 1
                    elif next_byte == b'\x16' and flag == 1:
                        data = serialPort.read(9)
                        Distance = struct.unpack('<f', data[5:9])[0]
                        if running:
                            self.current_data["Distance"] = Distance
                            self.queue.put({
                                "type": "distance",
                                "message": f"Distance: {Distance}\n",
                                "Distance": Distance
                            })
                            # Write the csv with all the data from the port
                            if self.csv_writer and all(k in self.current_data for k in ("timestamp", "HeartPhase", "BreathPhase", "TotalPhase", "Distance")):
                                current_time = datetime.datetime.strptime(self.current_data["timestamp"], "%H:%M:%S.%f")
                                if self.previous_timestamp is None:
                                    time_difference = 0.0  # First sample
                                else:
                                    time_difference = (current_time - self.previous_timestamp).total_seconds()
                                self.csv_writer.writerow([
                                    self.current_data["timestamp"],
                                    self.current_data["HeartPhase"],
                                    self.current_data["BreathPhase"],
                                    self.current_data["TotalPhase"],
                                    self.current_data["Distance"],
                                    time_difference
                                ])
                                self.previous_timestamp = current_time  # Update last timestamp
                                self.current_data = {}  # Clear buffer

            except Exception as e:
                self.queue.put({"type": "error", "message": f"Error reading serial data: {e}\n"})

    def update_ui(self):
        global running
        """Update the UI with new data from the queue."""
        
        while not self.queue.empty():
            item = self.queue.get()
            if item["type"] == "data":
                if running:
                    # Update labels
                    self.update_square_labels(item["HeartPhase"], item["BreathPhase"], item["TotalPhase"], None)
            elif item["type"] == "distance":
                if running:
                    # Update only the distance label
                    self.square_labels[3].config(text=f"Distance:\n{round(item['Distance'], decimal_round)}")

            # Log message (limit log size to 500 lines)
            self.text_log.config(state=tk.NORMAL)
            self.text_log.insert(tk.END, item["message"])
            if float(self.text_log.index('end-1c')) > 500:  # Keep max 500 lines
                self.text_log.delete("1.0", "2.0")
            self.text_log.config(state=tk.DISABLED)
            self.text_log.yview(tk.END)

        if running:
            # Update the graphs efficiently (redraw every 2 updates)
            self.update_counter += 1
            if self.update_counter % 2 == 0:
                self.line1.set_ydata(heart_phase_values)
                self.line2.set_ydata(breath_phase_values)
                self.canvas.draw_idle()  # Use draw_idle for better performance

        self.root.after(update_time, self.update_ui)

    def update_square_labels(self, HeartPhase, BreathPhase, TotalPhase, Distance):
        """Update the square labels with the latest values."""
        self.square_labels[0].config(text=f"Heart:\n{round(HeartPhase, decimal_round)}")
        self.square_labels[1].config(text=f"Breath:\n{round(BreathPhase, decimal_round)}")
        self.square_labels[2].config(text=f"Total:\n{round(TotalPhase, decimal_round)}")
        if Distance is not None:
            self.square_labels[3].config(text=f"Distance:\n{round(Distance, decimal_round)}")

    def on_close(self):
        """Handle application shutdown."""
        global running
        running = False
        if serialPort and serialPort.is_open:
            serialPort.close()
        if self.csv_file:
            self.csv_file.close()
        self.root.quit()
        self.root.destroy()

def create_gui():
    """Create and configure the GUI."""
    global root
    root = tk.Tk()
    root.geometry("1200x800")
    root.rowconfigure(0, weight=4)
    root.rowconfigure(1, weight=1)
    root.columnconfigure(0, weight=1)
    app = SerialReaderApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    return root

def main():
    """Main function to initialize and run the application."""
    initialize_serial()
    root = create_gui()
    root.mainloop()

if __name__ == "__main__":
    main()
