#include <QMdiSubWindow>
#include <cmath>
#include <cstddef>

#include "widget.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <math.h>
#include <qtimer.h>
#include <rtxi/rt.hpp>
#include <rtxi/rtos.hpp>

#include "ui_Membrane_Test_MainWindow.h"

membrane_test::Plugin::Plugin(Event::Manager* ev_manager)
    : Widgets::Plugin(ev_manager, std::string(membrane_test::MODULE_NAME))
{
  if (RT::OS::getFifo(this->Fifo, 100000) != 0) {
    ERROR_MSG("Unable to create Fifo for membrane_test plugin");
  }
}

membrane_test::Panel::Panel(QMainWindow* main_window,
                            Event::Manager* ev_manager)
    : Widgets::Panel(
          std::string(membrane_test::MODULE_NAME), main_window, ev_manager)
{
  setWhatsThis("Christini Lab Membrane Properties Probe");
  customizeGUI();
}

membrane_test::Component::Component(Widgets::Plugin* hplugin)
    : Widgets::Component(hplugin,
                         std::string(membrane_test::MODULE_NAME),
                         membrane_test::get_default_channels(),
                         membrane_test::get_default_vars())
    , fifo(dynamic_cast<membrane_test::Plugin*>(getHostPlugin())->getFifo())
{
  // The buffer will start with enough memory to handle 10000 values per
  // pulse width. In the event that the update frequency causes the
  // calculations to exceed this number we'll just have to bite the
  // bullet and potentially allocate additional memory when push_back is called.
  // This however should only happen once, and consecutive calls to execute
  // would not be a problem.... Unless you increase the frequency again.
  mp_data.reserve(10000);
}

void membrane_test::Component::execute()
{
  switch (this->getState()) {
    case RT::State::EXEC: {
      const int64_t current_time_ns =
          RT::OS::getTime() - this->measure_start_ns;
      const int period_section = (current_time_ns / (pulseWidth * 1000000)) % 2;
      switch (period_section) {
        case 0:
          writeoutput(0, (holdingVoltage + pulseAmp) * 1e-3);
          break;
        case 1:
          writeoutput(0, holdingVoltage * 1e-3);
          break;
        default:
          break;
      }
      const bool mp_collectData = acquire_data
          && (mp_stepsCount < mp_stepsTotal
              || mp_mode == mp_mode_t::CONTINUOUS);
      // Is this the first sample after a wave period? If so we should send data
      // and reset the current measurements to 0
      if (mp_collectData && (current_time_ns % pulseWidth) <= mp_period) {
        fifo->writeRT(mp_data.data(), sizeof(double) * mp_data.size());
        mp_data.clear();
      }
      if (mp_collectData) {
        // This try-catch block guarantees that our execution function does not
        // bring down the whole system because of unreasonable configuration.
        // It is better then crashing. Not to worry the data is preallocated and
        // it is unlikely that this would allocate. But if it does and fails
        // won't crash.
        try {
          mp_data.push_back(readinput(0));
        } catch (const std::bad_alloc& e) {
          ERROR_MSG(
              "membrane_test::Component::execute : Memory allocation failed "
              "after period change!");
          ERROR_MSG("Consider decreasing pulseWidth or increasing RT period");
          // Make sure we don't keep throwing.
          mp_data.clear();
          mp_stepsTotal = 0;
          acquire_data = false;
        }
      }
      break;
    }
    case RT::State::INIT:
      pulseAmp = getValue<double>(PULSE_AMP);
      holdingVoltage = getValue<double>(HOLDING_VOLTAGE);
      pulseWidth = getValue<int64_t>(PULSE_WIDTH);
      mp_stepsTotal = getValue<uint64_t>(TARGET_PULSE_COUNT);
      mp_stepsCount = 0;
      acquire_data = getValue<uint64_t>(ACQUIRE_ON) == 1;
      mp_mode = static_cast<mp_mode_t>(getValue<uint64_t>(MP_MODE));
      mp_period = RT::OS::getPeriod();
      setState(RT::State::PAUSE);
      break;
    case RT::State::MODIFY:
      pulseAmp = getValue<double>(PULSE_AMP);
      holdingVoltage = getValue<double>(HOLDING_VOLTAGE);
      pulseWidth = getValue<int64_t>(PULSE_WIDTH);
      std::cout << pulseWidth << std::endl;
      mp_stepsTotal = getValue<uint64_t>(TARGET_PULSE_COUNT);
      mp_stepsCount = 0;
      acquire_data = getValue<uint64_t>(ACQUIRE_ON) == 1;
      mp_mode = static_cast<mp_mode_t>(getValue<uint64_t>(MP_MODE));
      mp_period = RT::OS::getPeriod();
      setState(RT::State::EXEC);
      break;
    case RT::State::PERIOD:
      this->mp_period = RT::OS::getPeriod();
      setState(RT::State::EXEC);
      break;
    case RT::State::PAUSE:
      break;
    case RT::State::UNPAUSE:
      this->measure_start_ns = RT::OS::getTime();
      setState(RT::State::EXEC);
      break;
    default:
      break;
  }
}

void membrane_test::Panel::MP_Calculate()
{
  double Vpp = 0.0;
  if (mtUi.holdingVoltage1_button->isChecked()) {
    Vpp = mtUi.holdingVoltage1_spinBox->value();
  } else if (mtUi.holdingVoltage2_button->isChecked()) {
    Vpp = mtUi.holdingVoltage2_spinBox->value();
  } else {
    Vpp = mtUi.holdingVoltage3_spinBox->value();
  }

  // Taken from electrophys_plugin, written by Jonathan Bettencourt
  // In short, uses area under capacitive transient to calculate Cm by using
  // exponential curve fitting
  const size_t data_size = mp_data.size();
  uint64_t mp_stepsTotal = mtUi.mp_updatePeriod_spinBox->value();
  // Average has been taken before calling this function
  // for (size_t i = 0; i < data_size; ++i) {
  //  mp_data.at(i) /= mp_stepsTotal;
  //}

  double I1 = 0.0;
  for (auto i = static_cast<size_t>(round(data_size / 2 - ceil(data_size / 8)));
       i < data_size / 2;
       ++i)
  {
    I1 += mp_data.at(i);
  }
  I1 /= ceil(data_size / 8);

  double I2 = 0.0;
  for (auto i = static_cast<size_t>(round(data_size - ceil(data_size / 8)));
       i < data_size;
       ++i)
  {
    I2 += mp_data.at(i);
  }
  I2 /= ceil(data_size / 8);

  // Units?
  double dt = RT::OS::getPeriod() * 1e-6;

  double Q11 = NAN;
  double tau1 = NAN;
  {
    Q11 = 0.0;
    for (size_t i = 0; i < data_size / 2 - 1; ++i) {
      Q11 += dt * 1e-3 * (mp_data.at(i) + mp_data[i + 1] - 2 * I1) / 2;
    }
    Q11 = fabs(Q11);

    size_t xi = 0;
    for (; mp_data.at(xi) <= mp_data.at(xi + 1); ++xi) {
    }

    double sy = 0.0;
    double Y = mp_data.at(xi);
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data.at(xi) * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
    for (size_t i = xi + 1; i < data_size / 2; ++i) {
      sy += dt * 1e-3 * (mp_data.at(i - 1) + mp_data.at(i)) / 2;

      Y += mp_data.at(i);
      SY += sy;
      tSY += (i - xi) * dt * 1e-3 * sy;
      YSY += mp_data.at(i) * sy;
      SYSY += sy * sy;
      t += (i - xi) * dt * 1e-3;
      tt += ((i - xi) * dt * 1e-3) * ((i - xi) * dt * 1e-3);
      Yt += (i - xi) * dt * 1e-3 * mp_data.at(i);
    }

    std::array<double, 9> A = {
        static_cast<double>(data_size) / 2 - xi,
        SY,
        t,
        SY,
        SYSY,
        tSY,
        t,
        tSY,
        tt,
    };
    std::array<double, 3> B = {
        Y,
        YSY,
        Yt,
    };
    std::array<double, 9> V {};
    std::array<double, 3> S {};
    std::array<double, 3> x {};

    gsl_matrix_view a = gsl_matrix_view_array(A.data(), 3, 3);
    gsl_matrix_view b = gsl_matrix_view_array(V.data(), 3, 3);
    gsl_vector_view c = gsl_vector_view_array(S.data(), 3);
    gsl_vector_view d = gsl_vector_view_array(x.data(), 3);
    gsl_vector_view e = gsl_vector_view_array(B.data(), 3);

    gsl_linalg_SV_decomp(&a.matrix, &b.matrix, &c.vector, &d.vector);
    gsl_linalg_SV_solve(&a.matrix, &b.matrix, &c.vector, &e.vector, &d.vector);
    tau1 = fabs(1.0 / x[1]);
  }

  double Q12 = NAN;
  double tau2 = NAN;
  {
    Q12 = 0.0;
    for (size_t i = data_size / 2; i < data_size; ++i)
      Q12 += dt * 1e-3 * (mp_data.at(i) + mp_data[i + 1] - 2 * I2) / 2;
    Q12 = fabs(Q12);

    size_t xi = data_size / 2;
    for (; mp_data.at(xi) >= mp_data[xi + 1]; ++xi) {
      ;
    }

    double sy = 0.0;
    double Y = mp_data.at(xi);
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data.at(xi) * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
    for (size_t i = xi + 1; i < data_size; ++i) {
      sy += dt * 1e-3 * (mp_data[i - 1] + mp_data.at(i)) / 2;

      Y += mp_data.at(i);
      SY += sy;
      tSY += (i - xi) * dt * 1e-3 * sy;
      YSY += mp_data.at(i) * sy;
      SYSY += sy * sy;
      t += (i - xi) * dt * 1e-3;
      tt += ((i - xi) * dt * 1e-3) * ((i - xi) * dt * 1e-3);
      Yt += (i - xi) * dt * 1e-3 * mp_data.at(i);
    }

    std::array<double, 9> A = {
        static_cast<double>(data_size - xi),
        SY,
        t,
        SY,
        SYSY,
        tSY,
        t,
        tSY,
        tt,
    };
    std::array<double, 3> B = {
        Y,
        YSY,
        Yt,
    };
    std::array<double, 9> V {};
    std::array<double, 3> S {};
    std::array<double, 3> x {};

    gsl_matrix_view a = gsl_matrix_view_array(A.data(), 3, 3);
    gsl_matrix_view b = gsl_matrix_view_array(V.data(), 3, 3);
    gsl_vector_view c = gsl_vector_view_array(S.data(), 3);
    gsl_vector_view d = gsl_vector_view_array(x.data(), 3);
    gsl_vector_view e = gsl_vector_view_array(B.data(), 3);

    gsl_linalg_SV_decomp(&a.matrix, &b.matrix, &c.vector, &d.vector);
    gsl_linalg_SV_solve(&a.matrix, &b.matrix, &c.vector, &e.vector, &d.vector);
    tau2 = fabs(1.0 / x[1]);
  }

  double tau = (tau1 + tau2) / 2.0;
  double Q1 = (Q11 + Q12) / 2.0;
  double Q2 = fabs(I1 - I2) * tau;
  double Qt = Q1 + Q2;
  double rt = Vpp * 1e-3 / fabs(I1 - I2);

  ra = tau * Vpp * 1e-3 / Qt;
  rm = rt - ra;
  cm = Qt * rt / (Vpp * 1e-3 * rm);

  ra = round(ra * 1e-6 * 10) / 10;
  rm = round(rm * 1e-6 * 10) / 10;
  cm = round(cm * 1e12 * 10) / 10;
}

void membrane_test::Panel::customizeGUI()
{
  // Initialize Qt designer derived widget
  mtWindow = new QWidget(this);
  mtUi.setupUi(mtWindow);

  // Add newly created widget to layout of this widget

  auto* layout = new QVBoxLayout();
  setLayout(layout);
  layout->addWidget(mtWindow);

  // Set Ui refresh rate
  auto* timer = new QTimer(this);
  timer->start(100);  // 100ms refresh rate

  // Set timers
  rs_timer = new QTimer(this);
  rs_timer->setSingleShot(true);

  // Connect mtUi elements to slot functions
  // Resistance measurement
  QObject::connect(mtUi.pulse_button,
                   &QAbstractButton::toggled,
                   this,
                   &membrane_test::Panel::toggle_pulse);
  QObject::connect(
      mtUi.holdingVoltage1_button, SIGNAL(clicked()), this, SLOT(modify()));
  QObject::connect(
      mtUi.holdingVoltage2_button, SIGNAL(clicked()), this, SLOT(modify()));
  QObject::connect(
      mtUi.holdingVoltage3_button, SIGNAL(clicked()), this, SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage1_spinBox,
                   SIGNAL(valueChanged(int)),
                   this,
                   SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage2_spinBox,
                   SIGNAL(valueChanged(int)),
                   this,
                   SLOT(modify()));
  QObject::connect(mtUi.holdingVoltage3_spinBox,
                   SIGNAL(valueChanged(int)),
                   this,
                   SLOT(modify()));
  QObject::connect(
      mtUi.pulseAmp_spinBox, SIGNAL(valueChanged(int)), this, SLOT(modify()));
  QObject::connect(
      mtUi.pulseWidth_spinBox, SIGNAL(valueChanged(int)), this, SLOT(modify()));
  // Membrane properties
  QObject::connect(mtUi.mp_acquire_button,
                   SIGNAL(toggled(bool)),
                   this,
                   SLOT(toggle_mp_acquire(bool)));
  QObject::connect(mtUi.mp_updatePeriod_spinBox,
                   SIGNAL(valueChanged(int)),
                   this,
                   SLOT(modify()));
  QObject::connect(
      mtUi.mp_steps_spinBox, SIGNAL(valueChanged(int)), this, SLOT(modify()));
  QObject::connect(mtUi.mp_mode_comboBox,
                   SIGNAL(currentIndexChanged(int)),
                   this,
                   SLOT(modify()));
  // Acquire button is only enabled during pulse
  QObject::connect(mtUi.pulse_button,
                   SIGNAL(toggled(bool)),
                   mtUi.mp_acquire_button,
                   SLOT(setEnabled(bool)));
  // Update rate is only enabled when membrane property calculation is off
  QObject::connect(mtUi.mp_acquire_button,
                   SIGNAL(toggled(bool)),
                   mtUi.mp_updatePeriod_spinBox,
                   SLOT(setDisabled(bool)));
  // Display updates
  QObject::connect(
      timer, SIGNAL(timeout(void)), this, SLOT(update_rm_display(void)));
  QObject::connect(
      rs_timer, SIGNAL(timeout(void)), this, SLOT(resize_rm_text(void)));
  // QObject::connect(
  //     mp_timer, SIGNAL(timeout(void)), this, SLOT(update_mp_display(void)));

  resizeMe();
}

void membrane_test::Component::initialize()
{
  // Resistance measurement variables
  holdingVoltage = 0;
  pulseAmp = 10;
  pulseWidth = 20;

  // Membrane properties variables
  mp_stepsTotal = 500;
  mp_mode = mp_mode_t::SINGLE;
}

// Update parameter values based on Ui
void membrane_test::Panel::modify()
{
  Widgets::Plugin* hplugin = getHostPlugin();
  // Make sure real-time thread is not in the middle of execution
  getHostPlugin()->setComponentState(RT::State::PAUSE);

  // Resistance measurement
  if (mtUi.holdingVoltage1_button->isChecked()) {
    hplugin->setComponentParameter<double>(
        PARAMETER::HOLDING_VOLTAGE, mtUi.holdingVoltage1_spinBox->value());
  } else if (mtUi.holdingVoltage2_button->isChecked()) {
    hplugin->setComponentParameter<double>(
        PARAMETER::HOLDING_VOLTAGE, mtUi.holdingVoltage2_spinBox->value());
  } else {
    hplugin->setComponentParameter<double>(
        PARAMETER::HOLDING_VOLTAGE, mtUi.holdingVoltage3_spinBox->value());
  }

  hplugin->setComponentParameter<double>(PARAMETER::PULSE_AMP,
                                         mtUi.pulseAmp_spinBox->value());
  hplugin->setComponentParameter<int64_t>(PARAMETER::PULSE_WIDTH,
                                          mtUi.pulseWidth_spinBox->value());
  const uint64_t mp_stepsTotal = mtUi.mp_steps_spinBox->value();
  hplugin->setComponentParameter<uint64_t>(PARAMETER::TARGET_PULSE_COUNT,
                                           mp_stepsTotal);

  // Membrane properties
  const auto mp_mode =
      static_cast<mp_mode_t>(mtUi.mp_mode_comboBox->currentIndex());

  // Set update period minimum based on number of steps to be averaged

  min_updatePeriod = static_cast<int>(
      ceil(mp_stepsTotal * mtUi.pulseWidth_spinBox->value() / 1e3));
  mtUi.mp_updatePeriod_spinBox->setMinimum(std::max(1, min_updatePeriod));
  mp_updatePeriod = mtUi.mp_updatePeriod_spinBox->value();

  getHostPlugin()->setComponentState(RT::State::MODIFY);
}

// Toggle slot functions
void membrane_test::Panel::toggle_pulse(bool on)
{
  on ? rs_timer->start(mp_updatePeriod * 1000) : rs_timer->stop();
  getHostPlugin()->setComponentState(on ? RT::State::UNPAUSE
                                        : RT::State::PAUSE);
}

void membrane_test::Panel::toggle_mp_acquire(bool on)
{
  Widgets::Plugin* hplugin = getHostPlugin();
  const uint64_t acq = on ? 1UL : 0UL;
  hplugin->setComponentParameter(ACQUIRE_ON, acq);
}

// Update slot functions
// Resistance measurement value
void membrane_test::Panel::update_rm_display()
{
  if (getHostPlugin()->getComponentState() != RT::State::EXEC) {
    return;  // Return if pulse is not on
  }
  RT::OS::Fifo* fifo =
      dynamic_cast<membrane_test::Plugin*>(getHostPlugin())->getFifo();
  size_t num_bytes_read =
      fifo->read(mp_data.data(), sizeof(double) * mp_data.size());
  if (num_bytes_read == 0) {
    return;
  }
  if (mp_data.size() == num_bytes_read / sizeof(double)) {
    double value = 0.0;
    while (fifo->read(&value, sizeof(double)) > 0) {
      mp_data.push_back(value);
    }
  }

  // Calcualte the running average
  const int mp_stepsTotal = mtUi.mp_updatePeriod_spinBox->value();
  const size_t data_size = mp_data.size();
  ++pulse_count;
  if (mp_data_average.size() == 0 || pulse_count <= 1) {
    mp_data_average.resize(data_size);
    std::copy(mp_data.begin(), mp_data.end(), mp_data_average.begin());
  } else {
    for (auto indx = 0; indx < std::min(data_size, mp_data_average.size());
         ++indx)
    {
      // calculate average
      mp_data_average.at(indx) =
          (mp_data_average.at(indx) * (pulse_count - 1) + mp_data.at(indx))
          / pulse_count;
    }
  }

  // Calculate currents and dI
  double I_1 = 0.0;
  for (auto i = static_cast<size_t>(round(data_size / 2 - ceil(data_size / 8)));
       i < data_size / 2;
       ++i)
  {
    I_1 += mp_data.at(i);
  }
  I_1 /= ceil(data_size / 8);

  double I_2 = 0.0;
  for (auto i = static_cast<size_t>(round(data_size - ceil(data_size / 8)));
       i < data_size;
       ++i)
  {
    I_2 += mp_data.at(i);
  }
  I_2 /= ceil(data_size / 8);

  dI = (I_1 - I_2);

  const double pulseAmp = mtUi.pulseAmp_spinBox->value();
  double R = fabs((pulseAmp * 1e-3) / dI);  // Resistance calculation
  size_t exp = 0;  // Exponent of resistance

  if (R != INFINITY) {
    while (R >= 1e3) {
      R *= 1e-3;  // Reduce R by an order of magnitude
      exp++;  // Increase exponent counter
    }
  }

  QString RString;
  QString::asprintf("%7.5g", R);

  // Choose appropriate suffic based on exponent
  if (exp != 0) {
    if (exp == 1) {
      RString.append(" K").append(omega);
    } else if (exp == 2) {
      RString.append(" M").append(omega);
    } else if (exp == 3) {
      RString.append(" G").append(omega);

    } else {
      QString suffic;
      QString::asprintf(" * 1e%lu", 3 * exp);
    }
  } else {
    RString.append(" ").append(omega);
  }
  mtUi.resistance_valueLabel->setText(RString);
  if (pulse_count >= mtUi.mp_steps_spinBox->value()) {
    MP_Calculate();
    update_mp_display();
  }
}

void membrane_test::Panel::resize_rm_text()
{
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
  int fontSizeGuess = 28;  // Minimum font size
  for (;; ++fontSizeGuess) {
    testFont.setPointSize(fontSizeGuess);
    const QRect testRect = QFontMetrics(testFont).boundingRect(
        labelRect, Qt::AlignCenter, labelText);
    if (testRect.height() >= labelRect.height()
        || testRect.width() >= labelRect.width() || fontSizeGuess > 100)
    {
      break;
    }
  }

  labelFont.setPointSize(fontSizeGuess - 1);
  mtUi.resistance_valueLabel->setFont(labelFont);
}

// Membrane property values
void membrane_test::Panel::update_mp_display()
{
  mtUi.cm_valueLabel->setText(QString::number(cm).append(" pF"));
  mtUi.ra_valueLabel->setText(QString::number(ra).append(" M").append(omega));
  mtUi.rm_valueLabel->setText(QString::number(rm).append(" M").append(omega));

  if (mtUi.mp_mode_comboBox->currentIndex() == SINGLE) {
    mtUi.mp_acquire_button->setChecked(false);
  }
}

// Event handling
// TODO: Properly handle changes in period on the real-time side
void membrane_test::Plugin::receiveEvent(Event::Object* event)
{
  uint64_t cnt = 0;
  // When thread rate is changed, update rate dependent parameters
  switch (event->getType()) {
    case Event::RT_PERIOD_EVENT:
      // Number of loops for complete step (2x pulse width)
      // cnt = ((2.0 * pulseWidth) * 1e-3)
      //     / (RT::OS::getPeriod() * 1e-9);
      // Restart membrane properties data collection
      // if (mp_collectData)
      //   mp_collectData = false;
    default:
      break;
  }
}

///////// DO NOT MODIFY BELOW //////////
// The exception is if your plugin is not going to need real-time functionality.
// For this case just replace the craeteRTXIComponent return type to nullptr.
// RTXI will automatically handle that case and won't attach a component to the
// real time thread for your plugin.

std::unique_ptr<Widgets::Plugin> createRTXIPlugin(Event::Manager* ev_manager)
{
  return std::make_unique<membrane_test::Plugin>(ev_manager);
}

Widgets::Panel* createRTXIPanel(QMainWindow* main_window,
                                Event::Manager* ev_manager)
{
  return new membrane_test::Panel(main_window, ev_manager);
}

std::unique_ptr<Widgets::Component> createRTXIComponent(
    Widgets::Plugin* host_plugin)
{
  return std::make_unique<membrane_test::Component>(host_plugin);
}

Widgets::FactoryMethods fact;

extern "C"
{
Widgets::FactoryMethods* getFactories()
{
  fact.createPanel = &createRTXIPanel;
  fact.createComponent = &createRTXIComponent;
  fact.createPlugin = &createRTXIPlugin;
  return &fact;
}
};

//////////// END //////////////////////
