#include "Membrane_Test.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_linalg.h>

#include <main_window.h>

#include <QtGlobal>
#include <QtWidgets>

namespace {
class Membrane_Test_SyncEvent : public RT::Event {
 public:
  int callback() {
    return 0;
  }
};
}

// Create Module Instance
extern "C" Plugin::Object *createRTXIPlugin() {
  return new Membrane_Test::Module();
}

// Workspace
static Workspace::variable_t vars[] = {
  // Input
  { "Input Current (A)",
    "Input current (A) from target cell",
    Workspace::INPUT, },
  // Outputs
  { "Output Voltage (V)",
    "Output voltage (V) to target cell or internal input",
    Workspace::OUTPUT, },
};

// Number of variables in vars
static size_t num_vars = sizeof(vars) / sizeof(Workspace::variable_t);

Membrane_Test::Module::Module() :
    QWidget(MainWindow::getInstance()->centralWidget()),RT::Thread(0),
    Workspace::Instance("Membrane Test Module",
                        vars, num_vars) {

  // Build module GUI
  setWindowTitle(QString::number(getID()) +
                 " Membrane Test Module");
  createGUI();
  initialize(); // Initialize module variables
  modify();
  show();
}

Membrane_Test::Module::~Module() {
  // Make sure real-time thread is not in the middle of execution
  setActive(false);
  Membrane_Test_SyncEvent event;
  RT::System::getInstance()->postEvent(&event);
}

void Membrane_Test::Module::execute() {
  // First section (voltage step on)
  if (idx < (cnt / 2)) {
    // Only 2nd half of current used for resistance measurement
    if (idx >= (cnt / 4))
      I_1 += input(0); // Current during voltage on step

    // Voltage step on, convert from mV to V
    output(0) = (holdingVoltage + pulseAmp) * 1e-3;
  }
  else { // Second section (voltage step off)
    if (idx >= (3 * cnt) / 4)
      I_2 += input(0); // Current during voltage off step

    // Voltage step off, convert from mV to V
    output(0) = (holdingVoltage) * 1e-3;
  }

  // All data during voltage on step is used for membrane properties
  // calculation
  if (mp_collectData) {
    mp_data.at(idx) += input(0);
  }

  // Increment index, resetting to 0 when idx equals cnt
  if (!(++idx %= cnt)) {
    // Difference current between voltage on and off
    dI = (I_1 - I_2) / (cnt / 4);
    I_1 = I_2 = 0.0; // Reset current

    // If membrane properties calculation is on and calculation is not already
    // happening
    if (mp_on && !mp_dataFinished) {
      // Flag check ensures data collection starts at the beginning of a voltage
      // step
      if (!mp_collectData) {
        mp_stepsDone = 1;
        mp_collectData = true;
        // Reset membrane properties data
        mp_data.clear();
        mp_data.resize(cnt, 0);
      }
      // If the number of desired steps to average is reached
      else if (++mp_stepsDone > mp_stepsTotal) {
        mp_dataFinished = true;
      }
    }
  }

}

int Membrane_Test::Module::MP_Calculate(void) {
  double Vpp = pulseAmp;
  size_t data_size = cnt;
  if (data_size != mp_data.size()) // Check to make sure data size is correct
    return 1;

  // Taken from electrophys_plugin, written by Jonathan Bettencourt
  // In short, uses area under capacitive transient to calculate Cm by using
  // exponential curve fitting
  for (size_t i = 0; i < data_size; ++i)
    mp_data[i] /= mp_stepsTotal;

  double I1 = 0.0;
  for (size_t i = static_cast<size_t>(round(data_size /2 - ceil(data_size/8)));
       i < data_size / 2; ++i)
    I1 += mp_data[i];
  I1 /= ceil(data_size/8);

  double I2 = 0.0;
  for(size_t i = static_cast<size_t>(round(data_size - ceil(data_size / 8)));
      i < data_size;++i)
    I2 += mp_data[i];
  I2 /= ceil(data_size/8);

  double dt = RT::System::getInstance()->getPeriod() * 1e-6;

  double Q11;
  double tau1;
  {
    Q11 = 0.0;
    for(size_t i = 0; i < data_size / 2 - 1; ++i)
      Q11 += dt * 1e-3 * (mp_data[i] + mp_data[i + 1] - 2 * I1) / 2;
    Q11 = fabs(Q11);

    size_t xi = 0;
    for(; mp_data[xi] <= mp_data[xi + 1]; ++xi);

    double sy = 0.0;
    double Y = mp_data[xi];
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data[xi] * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
    for(size_t i = xi + 1; i < data_size / 2; ++i) {
      sy += dt * 1e-3 * (mp_data[i-1] + mp_data[i]) / 2;

      Y += mp_data[i];
      SY += sy;
      tSY += (i-xi) * dt * 1e-3 * sy;
      YSY += mp_data[i] * sy;
      SYSY += sy * sy;
      t += (i - xi) * dt * 1e-3;
      tt += ((i - xi) * dt * 1e-3) * ((i-xi) * dt * 1e-3);
      Yt += (i-xi) * dt * 1e-3 * mp_data[i];
    }

    double A[3 * 3] = {
      data_size / 2 - xi, SY, t,
      SY, SYSY, tSY,
      t, tSY, tt,
    };
    double B[3] = {
      Y,
      YSY,
      Yt,
    };
    double V[3 * 3];
    double S[3];
    double x[3];

    gsl_matrix_view a = gsl_matrix_view_array(A, 3, 3);
    gsl_matrix_view b = gsl_matrix_view_array(V, 3, 3);
    gsl_vector_view c = gsl_vector_view_array(S, 3);
    gsl_vector_view d = gsl_vector_view_array(x, 3);
    gsl_vector_view e = gsl_vector_view_array(B, 3);

    gsl_linalg_SV_decomp(
        &a.matrix,
        &b.matrix,
        &c.vector,
        &d.vector
                         );
    gsl_linalg_SV_solve(
        &a.matrix,
        &b.matrix,
        &c.vector,
        &e.vector,
        &d.vector
                        );
    tau1 = fabs(1.0 / x[1]);
  }

  double Q12;
  double tau2;
  {
    Q12 = 0.0;
    for(size_t i = data_size / 2; i < data_size; ++i)
      Q12 += dt * 1e-3 * (mp_data[i] + mp_data[i + 1] - 2 * I2) / 2;
    Q12 = fabs(Q12);

    size_t xi = data_size / 2;
    for(; mp_data[xi] >= mp_data[xi + 1]; ++xi);

    double sy = 0.0;
    double Y = mp_data[xi];
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data[xi] * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
    for(size_t i = xi + 1; i < data_size; ++i) {
      sy += dt * 1e-3 * (mp_data[i - 1] + mp_data[i]) / 2;

      Y += mp_data[i];
      SY += sy;
      tSY += (i - xi) * dt * 1e-3 * sy;
      YSY += mp_data[i] * sy;
      SYSY += sy * sy;
      t += (i-xi) * dt * 1e-3;
      tt += ((i-xi) * dt * 1e-3) * ((i-xi) * dt * 1e-3);
      Yt += (i-xi) * dt * 1e-3 * mp_data[i];
    }

    double A[3 * 3] = {
      data_size - xi, SY, t,
      SY, SYSY, tSY,
      t, tSY, tt,
    };
    double B[3] = {
      Y,
      YSY,
      Yt,
    };
    double V[3 * 3];
    double S[3];
    double x[3];

    gsl_matrix_view a = gsl_matrix_view_array(A, 3, 3);
    gsl_matrix_view b = gsl_matrix_view_array(V, 3, 3);
    gsl_vector_view c = gsl_vector_view_array(S,3);
    gsl_vector_view d = gsl_vector_view_array(x,3);
    gsl_vector_view e = gsl_vector_view_array(B,3);

    gsl_linalg_SV_decomp(
        &a.matrix,
        &b.matrix,
        &c.vector,
        &d.vector
                         );
    gsl_linalg_SV_solve(
        &a.matrix,
        &b.matrix,
        &c.vector,
        &e.vector,
        &d.vector
                        );
    tau2 = fabs(1.0 / x[1]);
  }

  double tau = (tau1 + tau2) / 2.0;

  double Q1 = (Q11 + Q12) / 2.0;
  double Q2 = fabs(I1 - I2) * tau;
  double Qt = Q1 + Q2;

  double Rt = Vpp * 1e-3 / fabs(I1 - I2);

  ra = tau * Vpp * 1e-3 / Qt;
  rm = rt - ra;
  cm = Qt * rt / (Vpp * 1e-3 * rm);

  ra = round(ra * 1e-6 * 10) / 10;
  rm = round(rm * 1e-6 * 10) / 10;
  cm = round(cm *1e12 * 10) / 10;

  return 0;
}

void Membrane_Test::Module::createGUI() {
  // Create subwindow and add it to main RTXI window
  subWindow = new QMdiSubWindow(MainWindow::getInstance());
  subWindow->setAttribute(Qt::WA_DeleteOnClose);
  subWindow->setWindowIcon(QIcon("/usr/local/lib/rtxi/RTXI-widget-icon.png"));
  subWindow->setWindowFlags(Qt::CustomizeWindowHint |
                            Qt::WindowCloseButtonHint |
                            Qt::WindowMinimizeButtonHint );
  MainWindow::getInstance()->createMdi(subWindow);
  // Set this widget to newly created subwindow
  subWindow->setWidget(this);

  // Initialize Qt designer derived widget
  mtWindow = new QWidget(this);
  mtUi.setupUi(mtWindow);

  // Add newly created widget to layout of this widget
  QVBoxLayout *layout = new QVBoxLayout(this);
  setLayout(layout);
  layout->addWidget(mtWindow);

  // Set Ui refresh rate
  QTimer *timer = new QTimer(this);
  timer->start(100); // 100ms refresh rate

  // Set timers
  rs_timer = new QTimer(this);
  rs_timer->setSingleShot(true);
  mp_timer = new QTimer(this);

  omega = QChar(0x3A9); // Greek letter omega for resistance

  // Connect mtUi elements to slot functions
  // Resistance measurement
  QObject::connect(mtUi.pulse_button, SIGNAL(toggled(bool)),
                   this, SLOT(toggle_pulse(bool)));
  QObject::connect(mtUi.holdingVoltage1_button, SIGNAL(clicked()),
                   this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage2_button, SIGNAL(clicked()),
                   this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage3_button, SIGNAL(clicked()),
                   this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage1_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage2_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage3_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.pulseAmp_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.pulseWidth_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  // Membrane properties
  QObject::connect(mtUi.mp_acquire_button, SIGNAL(toggled(bool)),
                   this, SLOT(toggle_mp_acquire(bool)));
  QObject::connect(mtUi.mp_updatePeriod_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.mp_steps_spinBox, SIGNAL(valueChanged(int)),
                   this, SLOT(modify()));
  QObject::connect(mtUi.mp_mode_comboBox, SIGNAL(currentIndexChanged(int)),
                   this, SLOT(modify()));
  // Acquire button is only enabled during pulse
  QObject::connect(mtUi.pulse_button, SIGNAL(toggled(bool)),
                   mtUi.mp_acquire_button, SLOT(setEnabled(bool)));
  // Update rate is only enabled when membrane property calculation is off
  QObject::connect(mtUi.mp_acquire_button, SIGNAL(toggled(bool)),
                   mtUi.mp_updatePeriod_spinBox, SLOT(setDisabled(bool)));
  // Display updates
  QObject::connect(timer, SIGNAL(timeout(void)),
                   this, SLOT(update_rm_display(void)));
  QObject::connect(rs_timer, SIGNAL(timeout(void)),
                   this, SLOT(resize_rm_text(void)));
  QObject::connect(mp_timer, SIGNAL(timeout(void)),
                   this, SLOT(update_mp_display(void)));

  subWindow->show();
  subWindow->adjustSize();
}

void Membrane_Test::Module::initialize() {
  // Resistance measurement variables
  holdingVoltage = 0;
  holdingVoltageOption_1 = 0;
  holdingVoltageOption_2 = -40;
  holdingVoltageOption_3 = -80;
  pulseAmp = 10;
  pulseWidth = 20;
  resistance = 0;

  // Number of loops for complete step (2x step width)
  cnt = ((2.0 * pulseWidth) * 1e-3) /
      (RT::System::getInstance()->getPeriod() * 1e-9);

  // Membrane properties variables
  mp_on = false;
  mp_dataFinished = false;
  mp_collectData = false;
  mp_updatePeriod = 30;
  mp_stepsTotal = 500;
  mp_stepsDone = 0;
  mp_mode = CONTINUOUS;
  cm = 0;
  ra = 0;
  rm = 0;
}

// Update parameter values based on Ui
void Membrane_Test::Module::modify() {
  bool active = getActive();
  // Make sure real-time thread is not in the middle of execution
  setActive(false);
  Membrane_Test_SyncEvent event;
  RT::System::getInstance()->postEvent(&event);

  // Restart data collection if parameters are changed
  mp_collectData = false;

  // Resistance measurement
  holdingVoltageOption_1 = mtUi.holdingVoltage1_spinBox->value();
  holdingVoltageOption_2 = mtUi.holdingVoltage2_spinBox->value();
  holdingVoltageOption_3 = mtUi.holdingVoltage3_spinBox->value();

  if (mtUi.holdingVoltage1_button->isChecked())
    holdingVoltage = holdingVoltageOption_1;
  else if (mtUi.holdingVoltage2_button->isChecked())
    holdingVoltage = holdingVoltageOption_2;
  else
    holdingVoltage = holdingVoltageOption_3;

  pulseAmp = mtUi.pulseAmp_spinBox->value();
  pulseWidth = mtUi.pulseWidth_spinBox->value();
  cnt = ((2.0 * pulseWidth) * 1e-3) /
      (RT::System::getInstance()->getPeriod() * 1e-9);

  // Membrane properties
  mp_stepsTotal = mtUi.mp_steps_spinBox->value();
  // Set update period minimum based on number of steps to be averaged
  mtUi.mp_updatePeriod_spinBox->setMinimum(mp_stepsTotal * pulseWidth / 1e3);
  mp_updatePeriod = mtUi.mp_updatePeriod_spinBox->value();
  mp_mode = static_cast<mp_mode_t>(mtUi.mp_mode_comboBox->currentIndex());

  setActive(active);
}

// Toggle slot functions
void Membrane_Test::Module::toggle_pulse(bool on) {
  setActive(on);
  if (!on) {
    if (mp_timer->isActive()) // Pause timer to prevent display updates
      mp_timer->stop();
    output(0) = 0;
  }
  else if (mp_on)
    // Start timer, convert s to ms
    mp_timer->start(mp_updatePeriod / 1e3);
}

void Membrane_Test::Module::toggle_mp_acquire(bool on) {
  if (on) {
    mp_collectData = false;
    mp_dataFinished = false;
    mp_on = true;
    // Start timer, convert s to ms
    mp_timer->start(mp_updatePeriod / 1e3);
  }
  else {
    if (mp_timer->isActive())
      mp_timer->stop();
    mp_on = false;
  }
}

// Update slot functions
// Resistance measurement value
void Membrane_Test::Module::update_rm_display() {
  if (!getActive())
    return; // Return if pulse is not on

  double R = fabs((pulseAmp * 1e-3) / dI); // Resistance calculation
  size_t exp = 0; // Exponent of resistance

  if (R != INFINITY)
    while (R >= 1e3) {
      R *= 1e-3; // Reduce R by an order of magnitude
      exp++; // Increase exponent counter
    }

  QString RString;
  RString.sprintf("%7.5g", R);

  // Choose appropriate suffic based on exponent
  if(exp) {
    if(exp == 1)
      RString.append(" K").append(omega);
    else if (exp == 2)
      RString.append(" M").append(omega);
    else if (exp == 3)
      RString.append(" G").append(omega);

    else {
      QString suffic;
      suffic.sprintf(" * 1e%lu", 3 * exp);
    }
  }
  else RString.append(" ").append(omega);

  mtUi.resistance_valueLabel->setText(RString);
}

void Membrane_Test::Module::resize_rm_text() {
  // Resize text if window was enlarged or shrunk
  // Grab current info of label
  QFont labelFont = mtUi.resistance_valueLabel->font();
  const QRect labelRect = mtUi.resistance_valueLabel->contentsRect();

  // Using placeholder text rather than actual text of label since boundingRect
  // cannot deal with RichText
  const QString labelText = "000.00 XXX";

  // Test increasing sizes until font size is too big, starting with a minimum
  // font size of 32
  QFont testFont(labelFont);
  int fontSizeGuess = 28; // Minimum font size
  for( ; ; ++fontSizeGuess) {
    testFont.setPointSize(fontSizeGuess);
    const QRect testRect =
        QFontMetrics(testFont).boundingRect(labelRect, Qt::AlignCenter,
                                            labelText);
    if (testRect.height() >= labelRect.height() ||
        testRect.width() >= labelRect.width() || fontSizeGuess > 100)
      break;
  }

  labelFont.setPointSize(fontSizeGuess - 1);
  mtUi.resistance_valueLabel->setFont(labelFont);
}

// Membrane property values
void Membrane_Test::Module::update_mp_display() {
  if (mp_dataFinished) { // Data has been collected
    int retval = MP_Calculate();

    if (retval) { // Error
      mtUi.cm_valueLabel->setText("Error");
      mtUi.ra_valueLabel->setText("Error");
      mtUi.rm_valueLabel->setText("Error");
    }
    else {
      mtUi.cm_valueLabel->setText(QString::number(cm).append(" pF"));
      mtUi.ra_valueLabel->
          setText(QString::number(ra).append(" M").append(omega));
      mtUi.rm_valueLabel->
          setText(QString::number(rm).append(" M").append(omega));
    }

    if (mp_mode == SINGLE) {
      mtUi.mp_acquire_button->setChecked(false);
      mp_timer->stop();
    }
    else { // Start next calculation
      // Reset flags
      mp_collectData = false;
      mp_dataFinished = false;
    }
  }
}

// Event handling
void Membrane_Test::Module::receiveEvent(const ::Event::Object *event) {
  // When thread rate is changed, update rate dependent parameters
  if(event->getName() == Event::RT_POSTPERIOD_EVENT) {
    // Number of loops for complete step (2x pulse width)
    cnt = ((2.0 * pulseWidth) * 1e-3) /
        (RT::System::getInstance()->getPeriod() * 1e-9);
    // Restart membrane properties data collection
    if (mp_collectData)
      mp_collectData = false;
  }
}

void Membrane_Test::Module::receiveEventRT(const ::Event::Object *event) {
  // When thread rate is changed, update rate dependent parameters
  if(event->getName() == Event::RT_POSTPERIOD_EVENT) {
    // Number of loops for complete step (2x pulse width)
    cnt = ((2.0 * pulseWidth) * 1e-3) /
        (RT::System::getInstance()->getPeriod() * 1e-9);
  }
}

// Resize resistance label text during window resizes
void Membrane_Test::Module::resizeEvent(QResizeEvent *) {
  rs_timer->start(50);
}

// Settings loading and saving
void Membrane_Test::Module::doLoad(const Settings::Object::State &s) {
  if (s.loadInteger("Maximized")) showMaximized();
  else if (s.loadInteger("Minimized")) showMinimized();

  if (s.loadInteger("W")) {
    subWindow->resize(s.loadInteger("W"), s.loadInteger("H"));
    parentWidget()->move(s.loadInteger("X"), s.loadInteger("Y"));
  }
}

void Membrane_Test::Module::doSave(Settings::Object::State &s) const {
  // Window settings
  if (subWindow->isMaximized())
    s.saveInteger("Maximized", 1);
  else if (subWindow->isMinimized())
    s.saveInteger("Minimized", 1);

  QPoint pos = subWindow->pos();
  s.saveInteger("X", pos.x());
  s.saveInteger("Y", pos.y());
  s.saveInteger("W", subWindow->width());
  s.saveInteger("H", subWindow->height());
}
