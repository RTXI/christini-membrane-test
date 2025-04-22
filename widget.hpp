
#include <string>

#include <rtxi/fifo.hpp>
#include <rtxi/widgets.hpp>

#include "ui_Membrane_Test_MainWindow.h"

// This is an generated header file. You may change the namespace, but
// make sure to do the same in implementation (.cpp) file
namespace membrane_test
{

enum mp_mode_t
{
  SINGLE = 0,
  CONTINUOUS = 1
};

constexpr std::string_view MODULE_NAME = "membrane_test";

enum PARAMETER : Widgets::Variable::Id
{
  // set parameter ids here
  HOLDING_VOLTAGE = 0,
  PULSE_AMP,
  PULSE_WIDTH,
  TARGET_PULSE_COUNT,
  ACQUIRE_ON,
  MP_MODE
};

inline std::vector<Widgets::Variable::Info> get_default_vars()
{
  return {
      {PARAMETER::HOLDING_VOLTAGE,
       "Holding Voltage",
       "Voltage (mV) level to hold for",
       Widgets::Variable::INT_PARAMETER,
       0.0},
      {PARAMETER::PULSE_AMP,
       "Pluse Amplitude",
       "The apmplitude of the applied pulse",
       Widgets::Variable::DOUBLE_PARAMETER,
       10.0},
      {PARAMETER::PULSE_WIDTH,
       "Pulse Width",
       "Width (ms) of the applied pulse",
       Widgets::Variable::INT_PARAMETER,
       20LL},
      {PARAMETER::TARGET_PULSE_COUNT,
       "Target Pulse Count",
       "Number of pulses to average over for membrane property calculation",
       Widgets::Variable::UINT_PARAMETER,
       0UL},
      {PARAMETER::ACQUIRE_ON,
       "Acquire data state",
       "1 to acquire membrane properties, 0 to not acquire membrane properties",
       Widgets::Variable::UINT_PARAMETER,
       0UL},
      {PARAMETER::MP_MODE,
       "Membrane Properties Mode",
       "Mode for acquiring membrane properties. 0 for single "
       "acquisition and 1 for continuous acquisition",
       Widgets::Variable::UINT_PARAMETER,
       0UL}};
}

inline std::vector<IO::channel_t> get_default_channels()
{
  return {
      {"Input Current (A)", "Input Current (A) from target cell", IO::INPUT},
      {"Output Voltage (V)",
       "Output voltage (V) to target cell or internal input",
       IO::OUTPUT}};
}

class Panel : public Widgets::Panel
{
  Q_OBJECT
public:
  Panel(QMainWindow* main_window, Event::Manager* ev_manager);
  void customizeGUI();  // Build and connect Ui elements

signals:
  void update_state();

public slots:
  void update_rm_display();  // Refresh resistance measurement display
  void resize_rm_text();  // Resizes resistance measurement font
  void update_pulse_button(); // synchronizes panel with component state
  void update_mp_display();  // Refresh membrane property values
  void modify() override;  // Update parameters
  void toggle_pulse(bool);  // Called when pulse button is pressed
  void toggle_mp_acquire(bool);  // Called when acquire button is pressed

  // Calculate membrane property values
  void MP_Calculate();

private:
  // We read current data from the realtime thread using this fifo
  std::vector<double> mp_data;
  std::vector<double> mp_data_average;
  Ui::Membrane_Test_UI mtUi;
  QWidget* mtWindow = nullptr;
  QMdiSubWindow* subWindow = nullptr;
  QTimer* rs_timer = nullptr;  // Resize timer
  QTimer* mp_timer = nullptr; // update values timer
  QTimer* state_timer = nullptr; // Update panel state
  double cm = 0;
  double ra = 0;
  double rm = 0;
  double dI = 0;
  int mp_updatePeriod = 0;  // Period for calc of membrane properties (s)
  int min_updatePeriod = 0;  // Minmum period for membrane property calc (s)
  int pulse_count = 0;  // Number of pulses since last acquisition
  const QChar omega = QChar(0x3A9);  // Greek letter omega for resistance
};

class Component : public Widgets::Component
{
public:
  explicit Component(Widgets::Plugin* hplugin);
  void execute() override;
  void initialize();

private:
  std::vector<double> mp_data;
  // We use this pipe/schannel to send data from the real time
  // to the ui thread without incurring uneccessary latencies
  RT::OS::Fifo* fifo = nullptr;

  // Resistance measurement variables
  uint64_t mp_stepsCount = 0;
  uint64_t mp_stepsTotal = 0;
  int64_t pulseWidth = 0;  // Width of voltage pulse
  int64_t measure_start_ns = 0;
  double holdingVoltage = 0;
  double pulseAmp = 0;  // Amplitude of voltage pulse
  int cycle_count = 0;  // Number of cycles since pulse started
  mp_mode_t mp_mode = SINGLE;
  bool acquire_data = false;
};

class Plugin : public Widgets::Plugin
{
public:
  explicit Plugin(Event::Manager* ev_manager);
  RT::OS::Fifo* getFifo() { return Fifo.get(); }
  // custom implementation for receiving events
  void receiveEvent(Event::Object* event) override;
  int64_t getPeriod() const { return period_ns; }

private:
  std::unique_ptr<RT::OS::Fifo> Fifo;
  int64_t period_ns = 0;
};

}  // namespace membrane_test
