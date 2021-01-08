<<<<<<< HEAD
from  django.urls import path
from . import views
from django.conf.urls import url

urlpatterns = [
    path('', views.home, name='home'),
    path('inverse', views.inverse, name='inverse'),
    path('direct', views.direct, name='direct')
=======
from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name="home"),
    path('direct', views.direct, name="direct"),
    path('inverse', views.inverse, name="inverse"),
>>>>>>> 4379fda8acb039391c8ec2bfe2ad5282d000b36b
]