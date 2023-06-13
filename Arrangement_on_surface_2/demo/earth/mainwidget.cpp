
#include "mainwidget.h"

#include <QMouseEvent>

#include <cmath>
#include <iostream>
#include <string>


MainWidget::~MainWidget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}


void MainWidget::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
{
  switch (e->button())
  {
  case Qt::LeftButton:
    m_left_mouse_button_down = flag;
    break;

  case Qt::MiddleButton:
    m_middle_mouse_button_down = flag;
    break;
  }
}
void MainWidget::mousePressEvent(QMouseEvent *e)
{
  set_mouse_button_pressed_flag(e, true);
  m_last_mouse_pos = QVector2D(e->position());
}
void MainWidget::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());
  const auto diff = current_mouse_pos - m_last_mouse_pos;

  if (m_left_mouse_button_down)
  {  
    const float rotation_scale_factor = 0.1f;
    const float theta_around_x = rotation_scale_factor * diff.y();
    const float theta_around_y = rotation_scale_factor * diff.x();
    m_camera.rotate(theta_around_x, theta_around_y);
  }
  else if(m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void MainWidget::mouseReleaseEvent(QMouseEvent *e)
{
  set_mouse_button_pressed_flag(e, false);
}
void MainWidget::timerEvent(QTimerEvent *)
{
  update();
}



void MainWidget::initializeGL()
{
  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_program();

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}

void MainWidget::init_camera()
{
  m_camera.set_pos(0, 0, 10);
}
void MainWidget::init_geometry()
{
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 3;
  m_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
}
void MainWidget::init_shader_program()
{
  m_shader_program.init();

  const char* vs = "shaders/vertex.glsl";
  const char* gs = "shaders/geometry.glsl";
  const char* fs = "shaders/fragment.glsl";
  m_shader_program.add_shader_from_file(vs, GL_VERTEX_SHADER);
  //m_shader_program.add_shader_from_file(gs, GL_GEOMETRY_SHADER);
  m_shader_program.add_shader_from_file(fs, GL_FRAGMENT_SHADER);

  m_shader_program.link();
  m_shader_program.validate();
}



void MainWidget::resizeGL(int w, int h)
{
  // Reset projection
  qreal aspect = qreal(w) / qreal(h ? h : 1);
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;
  m_camera.perspective(fov, aspect, z_near, z_far);
}


void MainWidget::paintGL()
{
  QMatrix4x4 model;
  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto mvp = projection * view * model;

  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  {
    m_shader_program.use();
    m_shader_program.set_uniform("MVP", mvp);
    
    m_sphere->draw();

    m_shader_program.unuse();
  }
}