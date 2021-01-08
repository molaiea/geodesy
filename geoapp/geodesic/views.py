from django.shortcuts import render
from .forms import finalform
# Create your views here.
def home(request) :
    return render(request, 'home.html')


def inverse(request) :
    return render(request, 'inverse.html')

def direct(request) :
    print(request.method)
    if request.method == 'POST':
        submit= request.POST.get("submit")
        form = finalform(request.POST)
        latitude, longitude = 0, 0
        print(form.is_valid())
        print(form.errors)
        if form.is_valid():
            latitude= form.cleaned_data.get("latitude")
            longitude= form.cleaned_data.get("longitude")
            a= form.cleaned_data.get("params")
            print(a)
            print(latitude)
        return render(request, 'direct.html', {'form': form, 'latitude': latitude, 'longitude': longitude})
    else:
        form = finalform()
    return render(request, 'direct.html', {'form': form})
    #chahia tayiba hh