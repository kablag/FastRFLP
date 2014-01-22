from django.conf.urls import patterns, url

from fastrflp import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index')
)