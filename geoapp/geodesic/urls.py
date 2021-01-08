from  django.urls import path
from . import views
from django.conf.urls import url

urlpatterns = [
    path('', views.home, name='home'),
    path('inverse', views.inverse, name='inverse'),
    path('direct', views.direct, name='direct')
]