
#include <string>

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
  PULSE_WIDTH
};

inline std::vector<Widgets::Variable::Info> get_default_vars()
{
  return {{PARAMETER::HOLDING_VOLTAGE,
           "Holding Voltage",
           "Voltage level to hold for",
           Widgets::Variable::INT_PARAMETER,
           0.0},
          {PARAMETER::PULSE_AMP,
           "Pluse Amplifier",
           "",
           Widgets::Variable::DOUBLE_PARAMETER,
           1.0},
          {PARAMETER::PULSE_WIDTH,
           "Pulse Width",
           "",
           Widgets::Variable::UINT_PARAMETER,
           0.0}};
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

public slots:
  void update_rm_display();  // Refresh resistance measurement display
  void resize_rm_text();  // Resizes resistance measurement font
  void update_mp_display();  // Refresh membrane property values
  void modify() override;  // Update parameters
  void toggle_pulse(bool);  // Called when pulse button is pressed
  void toggle_mp_acquire(bool);  // Called when acquire button is pressed
private:
  QWidget* mtWindow;
  QMdiSubWindow* subWindow;
  Ui::Membrane_Test_UI mtUi;
  QTimer* mp_timer;  // Membrane properties timer
  QTimer* rs_timer;  // Resize timer
  QChar omega;  // Greek letter omega for resistance
};

class Component : public Widgets::Component
{
public:
  explicit Component(Widgets::Plugin* hplugin);
  void execute() override;
  void initialize();
  int MP_Calculate();  // Calculate membrane property values

private:
  // Resistance measurement variables
  int holdingVoltage;
  int holdingVoltageOption_1;
  int holdingVoltageOption_2;
  int holdingVoltageOption_3;
  int pulseAmp;  // Amplitude of voltage pulse
  int pulseWidth;  // Width of voltage pulse
  double I_1;  // Current during positive pulse
  double I_2;  // Current during negative pulse
  double dI;  // Difference current
  double resistance;  // Calculated membrane resistance

  bool mp_on;
  bool mp_collectData;
  bool mp_dataFinished;
  std::vector<double> mp_data;  // Vector holding current data
  int mp_updatePeriod;  // Period at which calculation occurs
  int mp_stepsTotal;  // Number of steps to be averaged
  int mp_stepsDone;
  double cm;
  double ra;
  double rm;
  int idx;
  int cnt;

  mp_mode_t mp_mode;

  // Additional functionality needed for RealTime computation is to be placed
  // here
};

class Plugin : public Widgets::Plugin
{
public:
  explicit Plugin(Event::Manager* ev_manager);
  void receiveEvent(Event::Object* event) override;
};

}  // namespace membrane_test
