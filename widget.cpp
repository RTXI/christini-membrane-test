#include <QMdiSubWindow>
#include <QTextStream>
#include <fmt/format.h>
#include <QTimer>
#include <cmath>
#include <cstddef>

#include "widget.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <math.h>
#include <rtxi/rt.hpp>
#include <rtxi/rtos.hpp>

#include "ui_Membrane_Test_MainWindow.h"

membrane_test::Plugin::Plugin(Event::Manager* ev_manager)
    : Widgets::Plugin(ev_manager, std::string(membrane_test::MODULE_NAME))
{
  if (RT::OS::getFifo(this->Fifo, 100000) != 0) {
    ERROR_MSG("Unable to create Fifo for membrane_test plugin");
  }
  Event::Object event(Event::Type::RT_GET_PERIOD_EVENT);
  ev_manager->postEvent(&event);
  this->period_ns = std::any_cast<int64_t>(event.getParam("period"));
}

void membrane_test::Plugin::receiveEvent(Event::Object* event)
{
  // This is the exact same implementation as the default receiveEvent provided
  // by the widgets base class.. except we wish to store the period in userspace
  // so that the plugin can quickly calculate membrane properties without having
  // to consult the realtime system for the current period every time.
  switch (event->getType()) {
    case Event::Type::RT_PERIOD_EVENT:
      this->period_ns = std::any_cast<int64_t>(event->getParam("period"));
      if (this->getPanel() == nullptr) {
        this->getPanel()->signal_state_change(RT::State::PERIOD);
      }
      break;
    default:
      break;
  }
}

membrane_test::Panel::Panel(QMainWindow* main_window,
                            Event::Manager* ev_manager)
    : Widgets::Panel(
          std::string(membrane_test::MODULE_NAME), main_window, ev_manager)
{
  setWhatsThis("Christini Lab Membrane Properties Probe");
  mp_data.resize(10000, {0.0, false});
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
    const int64_t current_time_ns = RT::OS::getTime() - this->measure_start_ns;
    const int period_section = (current_time_ns / (pulseWidth * 1000000)) % 2;
      // generate the square wave
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
      // Check if a pulseWidth cycle has occurred
      const int64_t count = current_time_ns / (pulseWidth * 1000000 * 2);
      if (count > cycle_count) {
        cycle_count = count;
        fifo->writeRT(mp_data.data(), sizeof(measurement) * mp_data.size());
        mp_data.clear();
      }
      // This try-catch block guarantees that our execution function does not
      // bring down the whole system because of unreasonable configuration.
      // It is better then crashing. Not to worry the data is preallocated and
      // it is unlikely that this would allocate. But if it does and fails
      // won't crash.
      try {
        const double input = readinput(0);
        mp_data.emplace_back(measurement{input, !static_cast<bool>(period_section)});
      } catch (const std::bad_alloc& e) {
        ERROR_MSG(
            "membrane_test::Component::execute : Memory allocation failed "
            "after period change!");
        ERROR_MSG("Consider decreasing pulseWidth or increasing RT period");
        // Make sure we don't keep throwing.
        mp_data.clear();
        acquire_data = false;
      }
      break;
    }
    case RT::State::INIT:
      readinput(0);
      pulseAmp = getValue<double>(PULSE_AMP);
      holdingVoltage = getValue<double>(HOLDING_VOLTAGE);
      pulseWidth = getValue<int64_t>(PULSE_WIDTH);
      mp_stepsTotal = getValue<uint64_t>(TARGET_PULSE_COUNT);
      mp_stepsCount = 0;
      cycle_count = 0;
      acquire_data = getValue<uint64_t>(ACQUIRE_ON) == 1;
      mp_mode = static_cast<mp_mode_t>(getValue<uint64_t>(MP_MODE));
      setState(RT::State::PAUSE);
      break;
    case RT::State::MODIFY:
      readinput(0);
      pulseAmp = getValue<double>(PULSE_AMP);
      holdingVoltage = getValue<double>(HOLDING_VOLTAGE);
      pulseWidth = getValue<int64_t>(PULSE_WIDTH);
      mp_stepsTotal = getValue<uint64_t>(TARGET_PULSE_COUNT);
      mp_stepsCount = 0;
      cycle_count = 0;
      acquire_data = getValue<uint64_t>(ACQUIRE_ON) == 1UL;
      mp_mode = static_cast<mp_mode_t>(getValue<uint64_t>(MP_MODE));
      setState(RT::State::UNPAUSE);
      break;
    case RT::State::PERIOD:
      // discard input values
      readinput(0);
      setState(RT::State::EXEC);
      break;
    case RT::State::PAUSE:
      readinput(0);
      writeoutput(0, 0);
      break;
    case RT::State::UNPAUSE:
      this->measure_start_ns = RT::OS::getTime();
      cycle_count = 0;
      mp_data.clear();
      setState(RT::State::EXEC);
      break;
    default:
      break;
  }
}

void membrane_test::Panel::MP_Calculate()
{
  const double Vpp = mtUi.pulseAmp_spinBox->value();

  // Taken from electrophys_plugin, written by Jonathan Bettencourt
  // In short, uses area under capacitive transient to calculate Cm by using
  // exponential curve fitting
  const size_t data_size = mp_data_average.size();
  uint64_t mp_stepsTotal = mtUi.mp_updatePeriod_spinBox->value();
  // Average has been taken before calling this function
  // for (size_t i = 0; i < data_size; ++i) {
  //  mp_data_average.at(i) /= mp_stepsTotal;
  //}

  // Compute I1 and I2 using explicit integer math to avoid integer
  // division/ceil pitfalls and unsigned underflow for small data_size.
  double I1 = 0.0;
  double I2 = 0.0;
  if (data_size == 0) {
    // No data: early return with zeros
    ra = 0.0;
    rm = 0.0;
    cm = 0.0;
    return;
  }

  const size_t n = data_size;
  const size_t segment = static_cast<size_t>(std::max(1.0, std::ceil(static_cast<double>(n) / 8.0)));
  const size_t mid = n / 2;

  // I1: average over the segment ending at mid (exclusive)
  size_t start1 = (mid > segment) ? (mid - segment) : 0;
  size_t end1 = mid;  // exclusive
  size_t count1 = (end1 > start1) ? (end1 - start1) : 0;
  if (count1 > 0) {
    for (size_t i = start1; i < end1; ++i) {
      I1 += mp_data_average.at(i);
    }
    I1 /= static_cast<double>(count1);
  }

  // I2: average over the final segment
  size_t start2 = (n > segment) ? (n - segment) : 0;
  size_t end2 = n;  // exclusive
  size_t count2 = (end2 > start2) ? (end2 - start2) : 0;
  if (count2 > 0) {
    for (size_t i = start2; i < end2; ++i) {
      I2 += mp_data_average.at(i);
    }
    I2 /= static_cast<double>(count2);
  }

  // Units seem to be milliseconds
  // double dt = RT::OS::getPeriod() * 1e-6;
  double dt =
      dynamic_cast<membrane_test::Plugin*>(this->getHostPlugin())->getPeriod()
      * 1e-6;

  double Q11 = NAN;
  double tau1 = NAN;
  bool no_curvature1 = false;
  {
    Q11 = 0.0;
    if (mid >= 2) {
      for (size_t i = 0; i + 1 < mid; ++i) {
        Q11 += dt * 1e-3 * (mp_data_average.at(i) + mp_data_average.at(i + 1) - 2 * I1) / 2;
      }
    }
    Q11 = fabs(Q11);

    // the max value SHOULD be the point where the curve starts falling at
    // the beginning of the square wave.
  // Search only in the first half for the maximum
  auto first_half_end = (mp_data_average.begin() + mid);
  auto iter = std::max_element(mp_data_average.begin(), first_half_end);
  long xi = std::distance(mp_data_average.begin(), iter);

    double sy = 0.0;
    double Y = mp_data_average.at(xi);
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data_average.at(xi) * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
  for (size_t i = static_cast<size_t>(xi) + 1; i < mid; ++i) {
      sy += dt * 1e-3 * (mp_data_average.at(i - 1) + mp_data_average.at(i)) / 2;

      Y += mp_data_average.at(i);
      SY += sy;
      tSY += (i - xi) * dt * 1e-3 * sy;
      YSY += mp_data_average.at(i) * sy;
      SYSY += sy * sy;
      t += (i - xi) * dt * 1e-3;
      tt += ((i - xi) * dt * 1e-3) * ((i - xi) * dt * 1e-3);
      Yt += (i - xi) * dt * 1e-3 * mp_data_average.at(i);
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
    double x1 = x.at(1);
    if (!std::isfinite(x1) || fabs(x1) < 1e-12) {
      // No measurable curvature in fit (division would be invalid)
      no_curvature1 = true;
      tau1 = NAN;
    } else {
      tau1 = fabs(1.0 / x1);
    }
  }

  double Q12 = NAN;
  double tau2 = NAN;
  bool no_curvature2 = false;
  {
    Q12 = 0.0;
    if (mid + 1 < n) {
      for (size_t i = mid; i + 1 < n; ++i) {
        Q12 += dt * 1e-3 * (mp_data_average.at(i) + mp_data_average.at(i + 1) - 2 * I2) / 2;
      }
    }
    Q12 = fabs(Q12);

    // the min value SHOULD be the point where the curve starts rising at
    // the middle of the square wave.
  auto second_half_begin = (mp_data_average.begin() + mid);
  auto iter = std::min_element(second_half_begin, mp_data_average.end());
  long xi = std::distance(mp_data_average.begin(), iter);

    double sy = 0.0;
    double Y = mp_data_average.at(xi);
    double SY = sy;
    double tSY = 0.0;
    double YSY = mp_data_average.at(xi) * sy;
    double SYSY = sy * sy;
    double t = 0.0;
    double tt = 0.0;
    double Yt = 0.0;
  for (size_t i = static_cast<size_t>(xi) + 1; i < n; ++i) {
      sy += dt * 1e-3 * (mp_data_average.at(i - 1) + mp_data_average.at(i)) / 2;

      Y += mp_data_average.at(i);
      SY += sy;
      tSY += (i - xi) * dt * 1e-3 * sy;
      YSY += mp_data_average.at(i) * sy;
      SYSY += sy * sy;
      t += (i - xi) * dt * 1e-3;
      tt += ((i - xi) * dt * 1e-3) * ((i - xi) * dt * 1e-3);
      Yt += (i - xi) * dt * 1e-3 * mp_data_average.at(i);
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
    double x1 = x.at(1);
    if (!std::isfinite(x1) || fabs(x1) < 1e-12) {
      // No measurable curvature in fit (division would be invalid)
      no_curvature2 = true;
      tau2 = NAN;
    } else {
      tau2 = fabs(1.0 / x1);
    }
  }

  double tau = NAN;
  if (std::isfinite(tau1) && std::isfinite(tau2)) {
    tau = (tau1 + tau2) / 2.0;
  }
  double Q1 = (Q11 + Q12) / 2.0;
  double deltaI = fabs(I1 - I2);
  double Q2 = deltaI * tau;
  double Qt = Q1 + Q2;

  double rt = NAN;
  if (deltaI > 1e-12) {
    rt = Vpp * 1e-3 / deltaI;
  }

  // Protect against division by zero
  if (Qt <= 0 || !std::isfinite(Qt) || !std::isfinite(tau) || !std::isfinite(rt)) {
    ra = 0.0;
    rm = 0.0;
    cm = 0.0;
    // If the cause was lack of curvature in the exponential fit, warn the
    // user so they know the measurement couldn't produce reliable values.
    if (no_curvature1 || no_curvature2) {
      ERROR_MSG("MP_Calculate: Fit reported no measurable curvature; ra/rm/cm set to 0. Check data quality or increase fitting window.");
    }
  } else {
    ra = tau * Vpp * 1e-3 / Qt;
    rm = rt - ra;
    if (rm == 0.0) {
      cm = 0.0;
    } else {
      cm = Qt * rt / (Vpp * 1e-3 * rm);
    }
  }

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

  // Set timers
  rs_timer = new QTimer(this);
  mp_timer = new QTimer(this);

  // Connect mtUi elements to slot functions
  // Resistance measurement
  QObject::connect(mtUi.pulse_button,
                   &QAbstractButton::toggled,
                   this,
                   &membrane_test::Panel::toggle_pulse);
  QObject::connect(mtUi.holdingVoltage1_button,
                   &QRadioButton::clicked,
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.holdingVoltage2_button,
                   &QRadioButton::clicked,
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.holdingVoltage3_button,
                   &QRadioButton::clicked,
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.holdingVoltage1_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.holdingVoltage2_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.holdingVoltage3_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.pulseAmp_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.pulseWidth_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  // Membrane properties
  QObject::connect(mtUi.mp_acquire_button,
                   &QPushButton::toggled,
                   this,
                   &membrane_test::Panel::toggle_mp_acquire);
  QObject::connect(mtUi.mp_updatePeriod_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.mp_steps_spinBox,
                   QOverload<int>::of(&QSpinBox::valueChanged),
                   this,
                   &membrane_test::Panel::modify);
  QObject::connect(mtUi.mp_mode_comboBox,
                   QOverload<int>::of(&QComboBox::currentIndexChanged),
                   this,
                   &membrane_test::Panel::modify);
  // Acquire button is only enabled during pulse
  QObject::connect(mtUi.pulse_button,
                   &QPushButton::toggled,
                   mtUi.mp_acquire_button,
                   &membrane_test::Panel::setEnabled);
  // Update rate is only enabled when membrane property calculation is off
  QObject::connect(mtUi.mp_acquire_button,
                   &QPushButton::toggled,
                   mtUi.mp_updatePeriod_spinBox,
                   &membrane_test::Panel::setDisabled);
  // Display updates
  QObject::connect(rs_timer,
                   &QTimer::timeout,
                   this,
                   &membrane_test::Panel::update_rm_display);
  QObject::connect(
      rs_timer, &QTimer::timeout, this, &membrane_test::Panel::resize_rm_text);
  QObject::connect(mp_timer,
                   &QTimer::timeout,
                   this,
                   &membrane_test::Panel::update_pulse_button);
  mp_timer->start(1000);
  // rs_timer->start(100);
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
  const RT::State::state_t prev_state = hplugin->getComponentState();
  // Make sure real-time thread is not in the middle of execution
  hplugin->setComponentState(RT::State::PAUSE);

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

  hplugin->setComponentParameter(ACQUIRE_ON,
                                 mtUi.mp_acquire_button->isDown() ? 1UL : 0UL);
  // Membrane properties
  const auto mp_mode =
      static_cast<mp_mode_t>(mtUi.mp_mode_comboBox->currentIndex());

  // Set update period minimum based on number of steps to be averaged

  min_updatePeriod = static_cast<int>(
      ceil(mp_stepsTotal * mtUi.pulseWidth_spinBox->value() / 1e3));
  mtUi.mp_updatePeriod_spinBox->setMinimum(std::max(1, min_updatePeriod));
  mp_updatePeriod = mtUi.mp_updatePeriod_spinBox->value();

  hplugin->setComponentState(RT::State::MODIFY);
  hplugin->setComponentState(prev_state);
}

// Toggle slot functions
void membrane_test::Panel::toggle_pulse(bool on)
{
  auto* hplugin = this->getHostPlugin();
  on ? rs_timer->start(mtUi.pulseWidth_spinBox->value()) : rs_timer->stop();
  getHostPlugin()->setComponentState(on ? RT::State::UNPAUSE
                                        : RT::State::PAUSE);
}

void membrane_test::Panel::toggle_mp_acquire(bool on)
{
  Widgets::Plugin* hplugin = getHostPlugin();
  pulse_count = 0;
  const uint64_t acq = on ? 1UL : 0UL;
  hplugin->setComponentParameter(ACQUIRE_ON, acq);
  hplugin->setComponentState(RT::State::MODIFY);
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
  int64_t num_bytes_read =
      fifo->read(mp_data.data(), sizeof(measurement) * mp_data.size());
  if (num_bytes_read <= 0) {
    return;
  }
  if (mp_data.size() == num_bytes_read / sizeof(measurement)) {
    measurement value = {0.0, false};
    while (fifo->read(&value, sizeof(measurement)) > 0) {
      mp_data.push_back(value);
      num_bytes_read += sizeof(double);
    }
  } else {
    mp_data.resize(num_bytes_read / sizeof(measurement));
  }

  const size_t data_size = num_bytes_read / sizeof(measurement);
  // Calculate currents and dI
  double I_1 = 0.0;
  int first_count = 0;
  double I_2 = 0.0;
  int second_count = 0;
  for (const auto& value : mp_data) {
    value.on ? I_1 += value.current : I_2 += value.current;
    value.on ? ++first_count : ++second_count;
  }
  I_1 /= first_count;
  I_2 /= second_count;

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

  // We want the decimal point to remain in the same column regardless of
  // value or unit. Use a fixed-width numeric field: 3 integer digits,
  // a decimal point, and 2 fractional digits -> total width = 6.
  // Pad with zeros when necessary so the visual column of the '.' remains
  // constant. Use fmt for consistent fixed formatting (no exponential).
  const int int_width = 3;
  const int precision = 2;  // at most 2 decimal places
  const int total_width = int_width + 1 + precision;  // e.g. 3 + '.' + 2 = 6

  const std::string formatted = fmt::format(
      "{:0{w}.{p}f}", R, fmt::arg("w", total_width), fmt::arg("p", precision));

  // Build unit string without leading space so we can layout numeric and unit
  // independently in HTML and guarantee the decimal point stays fixed.
  QString unit;
  if (exp == 0) {
    unit = omega;  // plain ohms (just the omega symbol)
  } else if (exp == 1) {
    unit = QString("K").append(omega);
  } else if (exp == 2) {
    unit = QString("M").append(omega);
  } else if (exp == 3) {
    unit = QString("G").append(omega);
  } else {
    unit = QString::fromStdString(fmt::format("* 1e{}", 3 * exp)).append(omega);
  }

  // Use the label's font metrics to compute a fixed pixel width for the
  // numeric field (based on the maximum pattern "000.00") so the decimal
  // point remains vertically aligned regardless of unit length. Render as
  // simple HTML with two spans: numeric (fixed width, right-aligned) and
  // unit (left of it). QLabel supports a subset of HTML/CSS which works
  // for this use-case.
  QFont labelFont = mtUi.resistance_valueLabel->font();
  QFontMetrics fm(labelFont);
  const QString maxPattern = QString::fromStdString(fmt::format("{:0{w}.{p}f}", 0.0,
                                                              fmt::arg("w", total_width),
                                                              fmt::arg("p", precision)));
  const int numWidthPx = fm.horizontalAdvance(maxPattern);

  const QString html = QString("<span style='display:inline-block; width:%1px; text-align:right;'>%2</span><span style='margin-left:6px;'>%3</span>")
                           .arg(numWidthPx)
                           .arg(QString::fromStdString(formatted))
                           .arg(unit);

  mtUi.resistance_valueLabel->setText(html);

  // Calcualte the running average
  const int mp_stepsTotal = mtUi.mp_updatePeriod_spinBox->value();
  bool last_on = false;
  size_t indx = 0;
  if (mp_data_average.size() == 0 || pulse_count <= 1) {
    mp_data_average.reserve(data_size);
    for (size_t i = 0; i < data_size; ++i) {
      if (mp_data.at(i).on && !last_on) {
        indx = 0;
        pulse_count++;
      }
      if (indx < mp_data_average.size()) {
        mp_data_average.at(indx) += mp_data.at(i).current;
      } else {
        mp_data_average.push_back(mp_data.at(i).current);
      }
      last_on = mp_data.at(i).on;
      indx++;
    }
  } else {
    for (auto i = 0; i < data_size; ++i) {
      if (mp_data.at(i).on && !last_on) {
        indx = 0;
        pulse_count++;
      }
      // calculate average
      if (mp_data_average.size() <= indx) {
        mp_data_average.push_back(mp_data.at(i).current);
      } else {
        mp_data_average.at(indx) = (mp_data_average.at(indx) * (pulse_count - 1)
                                    + mp_data.at(i).current)
            / pulse_count;
      }
    }
  }

  if (this->mtUi.mp_acquire_button->isChecked()
      && pulse_count >= mtUi.mp_steps_spinBox->value())
  {
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

void membrane_test::Panel::update_pulse_button()
{
  if (getHostPlugin() == nullptr) {
    return;
  }
  // Make sure real-time thread is not in the middle of execution
  Widgets::Plugin* hplugin = getHostPlugin();
  const RT::State::state_t state = getHostPlugin()->getComponentState();
  if (state == RT::State::UNDEFINED) {
    return;
  }
  const bool paused = state != RT::State::PAUSE;
  mtUi.pulse_button->setChecked(paused);
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
