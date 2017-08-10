#ifndef Membrane_Test_H
#define Membrane_Test_H

#include "Membrane_Test_MainWindow_UI.h"

#include <rt.h>
#include <settings.h>
#include <workspace.h>
#include <event.h>
#include <plugin.h>

#include <QtGlobal>
#include <QtWidgets>

namespace Membrane_Test {
class Module: public QWidget, public RT::Thread, public Plugin::Object,
              public Workspace::Instance, public Event::Handler,
              public Event::RTHandler {
  Q_OBJECT // Required macro for QT slots

 public:
  Module();
  ~Module();
  void execute(); // Function run at every RTXI loop
  void receiveEvent(const ::Event::Object *); // Receive non-RT event
  void receiveEventRT(const ::Event::Object *); // Receive RT event

 public slots:
  void update_rm_display(); // Refresh resistance measurement display
  void resize_rm_text(); // Resizes resistance measurement font
  void update_mp_display(); // Refresh membrane property values
  void modify(); // Update parameters
  void toggle_pulse(bool); // Called when pulse button is pressed
  void toggle_mp_acquire(bool); // Called when acquire button is pressed

 private:
  // Ui elements
  QWidget *mtWindow;
  QMdiSubWindow *subWindow;
  Ui::Membrane_Test_UI mtUi;
  QTimer *mp_timer; // Membrane properties timer
  QTimer *rs_timer; // Resize timer
  QChar omega; // Greek letter omega for resistance

  // Counters
  int idx;
  int cnt;

  // Resistance measurement variables
  int holdingVoltage;
  int holdingVoltageOption_1;
  int holdingVoltageOption_2;
  int holdingVoltageOption_3;
  int pulseAmp; // Amplitude of voltage pulse
  int pulseWidth; // Width of voltage pulse
  double I_1; // Current during positive pulse
  double I_2; // Current during negative pulse
  double dI; // Difference current
  double resistance; // Calculated membrane resistance

  // Membrane properties
  bool mp_on;
  bool mp_collectData;
  bool mp_dataFinished;
  std::vector<double> mp_data; // Ve
  double mp_updateRate; // Rate (Hz) at which calculation occurs
  int mp_stepsTotal; // Number of steps to be averaged
  int mp_stepsDone;
  double rt;
  double cm;
  double ra;
  double rm;

  enum mp_mode_t { SINGLE = 0, CONTINUOUS = 1 };
  mp_mode_t mp_mode;

  // Module functions
  void createGUI(); // Build and connect Ui elements
  void initialize(); // Initialize module parameters
  int MP_Calculate(); // Calculate membrane property values

 protected:
  void doLoad(const Settings::Object::State &);
  void doSave(Settings::Object::State &) const;
  // Reimplement resize event function for resizing of resistance label text
  virtual void resizeEvent(QResizeEvent *);
}; // Class Module
}; // Namespace Membrane_Test

#endif // Membrane_Test_H
