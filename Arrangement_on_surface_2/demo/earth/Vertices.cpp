
#include "Vertices.h"


Vertices::Vertices(const std::vector<QVector3D>& vertices)
{
  initializeOpenGLFunctions();

  auto& vertex_data = vertices;
  m_num_indices = vertices.size();

  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);

  // Vertex Buffer
  glGenBuffers(1, &m_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  auto vertex_buffer_size = sizeof(QVector3D) * vertex_data.size();
  auto vertex_buffer_data = reinterpret_cast<const void*>(vertex_data.data());
  glBufferData(GL_ARRAY_BUFFER,
               vertex_buffer_size,
               vertex_buffer_data,
               GL_STATIC_DRAW);

  // Position Vertex-Attribute
  GLint position_attrib_index = 0;
  const void* position_offset = 0;
  GLsizei stride = 0;
  glVertexAttribPointer(position_attrib_index,
                        3,
                        GL_FLOAT, GL_FALSE,
                        stride,
                        position_offset);
  glEnableVertexAttribArray(position_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void Vertices::draw()
{
  glBindVertexArray(m_vao);
  {
    glDrawArrays(GL_POINTS, 0, m_num_indices);
  }
  glBindVertexArray(0);
}
